%%% change line 7 p_list for range
%%% change line 44 for E/I balance or only E,
%%% change line 14 and 15 for steps

clearvars;
p_list = logspace(-1.3,4,20);     % balance
%p_list = logspace(-1.3,2,35);    % only Excitation
repeat = 3;
IMSE_raw = NaN*zeros(length(p_list),repeat);
SNR_indirect_raw = NaN*zeros(length(p_list),repeat);
SNR_indirect_baseline_raw = NaN*zeros(length(p_list),repeat);
SNR_indirect_Lesica_raw = NaN*zeros(length(p_list),repeat);
SNR_KL_full_raw = NaN*zeros(length(p_list),repeat);
SNR_KL_woH_raw = NaN*zeros(length(p_list),repeat);

tot_t = 1e6;
bin = 1;   %ms
ddt = bin;
V_E = 0.02;
V_I = V_E;
V_th = 1;
V_reset = 0;
tau_E = 1e-3;       % 1ms
tau_I = tau_E;
tau_M = 20;
dt = 1;

maxSig = 0.9e-1;
signalType = 3; % 1 for no signal, 2 for square wave, 3 for gamma
I_per = zeros(1,1e3);
if signalType == 2
    I_per(1:1000) = 1;
    I_per = I_per/max(I_per)*maxSig;
end
if signalType == 3
    pd = makedist('InverseGaussian','mu',2e2,'lambda',6e2);
    I_per(201:1000) = pdf(pd,[1:800]);
    I_per = I_per/max(I_per)*maxSig;
end
%plot(I_per);
repnum = ceil(tot_t/length(I_per));
I = repmat(I_per,1,repnum);
I = I(1:tot_t);
I_noInp = zeros(1,tot_t);

for ii = 1:length(p_list)
    p = p_list(ii);
    q = 0;
    q = p;
    for jj = 1:repeat
        [ii,jj]
        % Simulation
        [ISI,spike_timing,y_sparse,V,inputE,inputI] = ... 
            GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt);
        y = full(y_sparse);

        T = 1e3;
        N = 100;
        tot_N = (tot_t/length(I_per));
        ST = (reshape(y,[],tot_N)');
        fr = sum(ST);
        if ddt ~= 1
            fr2 = reshape(fr,ddt,[]);
            fr2 = sum(fr2);
        else
            fr2 = fr;
        end

        % MISE
        target = fr2/tot_N/ddt;
        groundTrue = I_per;
        %baseline = mean(target(1:200));
        %fun = @(k) sum(((target-baseline)*k - groundTrue).^2)
        %x0 = 1;
        fun = @(k) sum(((target-k(2))*k(1) - groundTrue).^2)
        x0 = [1,0.2];
        [x,fval] = fminunc(fun,x0);
        IMSE_raw(ii,jj) = fval;
        
        SNR_indirect_raw(ii,jj) = var(mean(ST)) / mean(var(ST));
        SNR_indirect_baseline_raw(ii,jj) = var(mean(ST(:,201:800))) / mean(var(ST(:,[1:200,801:1000])));
        SNR_indirect_Lesica_raw(ii,jj) = var(mean(ST)) / mean( var( (ST-mean(ST))' )' );
        

        %%% GLM with stimulus filter
        y_glm = reshape(ST',1,[]);
        y_glm(find(y_glm>1))=1;
        I_AllInputs = I;
        plotFlag = 0;
        plotKS = 0;
        dt = 1;
        fit_k = 1;

        nkt = 100; % number of ms in stim filter
        kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
        kbasprs.ncos = 5; % number of raised-cosine vectors to use
        kbasprs.kpeaks = [1 round(nkt/1.5)];  % position of first and last bump (relative to identity bumps)
        kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
        %%% basis functions for post-spike kernel
        ihbasprs.ncols = 5;  % number of basis vectors for post-spike kernel
        hPeaksMax = 200;
        ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
        ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
        ihbasprs.absref = 0; % absolute refractory period, in ms

        [k, h, dc, prs, kbasis, hbasis, stats] = ...
            fit_glm(I_AllInputs',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
        [pvalue, rate, h_output, k_output] = KStest(y_glm, h', I_AllInputs, k', dc, plotKS);
        nlogl_0HS = -sum( log(rate).*y_glm - rate );    % 0HS : baseline, history, stimulus ALL included
        %nlogl_satured = sum( y_glm );
        nlogl_satured = 0;
        Dev_0HS = 2* (nlogl_0HS - nlogl_satured);

        % without stimulus filter
        fit_k = 0;

        [k, h, dc, prs, kbasis, hbasis, stats] = ...
            fit_glm(zeros(1,length(I_AllInputs))',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
        k = zeros(1,length(k))';
        [pvalue, rate, h_output, k_output] = KStest(y_glm, h', I_AllInputs, k', dc, plotKS);
        nlogl_0H = -sum( log(rate).*y_glm - rate );    % 0H : baseline, history included
        Dev_0H = 2* (nlogl_0H - nlogl_satured);
        % without history filter
        fit_k = 2;

        [k, h, dc, prs, kbasis, hbasis, stats] = ...
            fit_glm(I_AllInputs',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
        [pvalue, rate, h_output, k_output] = KStest(y_glm, h', I_AllInputs, k', dc, plotKS);
        nlogl_0S = -sum( log(rate).*y_glm - rate );    % 0S : baseline, stimulus included
        Dev_0S = 2* (nlogl_0S - nlogl_satured);

        % without both history filter and stimulus filter, only baseline
        rate0 = ones(1,length(y_glm)) * sum(y_glm)/length(y_glm);
        nlogl_0 = -sum( log(rate0).*y_glm - rate0 );    % 0 : baseline included
        Dev_0 = 2* (nlogl_0 - nlogl_satured);

        % 
        SNR_glmS = (Dev_0H - Dev_0HS +kbasprs.ncos )...
            /(Dev_0HS +1+kbasprs.ncos+ihbasprs.ncols );
        SNR_glmH = (Dev_0S - Dev_0HS +ihbasprs.ncols )...
            /(Dev_0HS +1+kbasprs.ncos+ihbasprs.ncols );
        SNR_glmS2 = (Dev_0 - Dev_0S +kbasprs.ncos )...
            /(Dev_0S +1+kbasprs.ncos );
        db_glmS = 10*log10(SNR_glmS);
        db_glmS2 = 10*log10(SNR_glmS2);
        db_glmH = 10*log10(SNR_glmH);

        EPEsignal = sum( (rate0 - rate).^2 );
        EPEnoise = sum( (y_glm - rate).^2 ) ;
        EPEsignal2 = sum( (y_glm - rate0).^2 ) - EPEnoise;
        SNR_EPE = EPEsignal/EPEnoise;
        db_EPE = 10*log10(SNR_EPE);


        SNR_KL_full_raw(ii,jj) = SNR_glmS;
        SNR_KL_woH_raw(ii,jj) = SNR_glmS2;
    end
end
%%
figure;
hold on
errorbar(p_list*1e3/50,mean(IMSE_raw,2),std(IMSE_raw')');
errorbar(p_list*1e3/50,mean(SNR_indirect_raw,2),std(SNR_indirect_raw')');
errorbar(p_list*1e3/50,mean(SNR_indirect_baseline_raw,2),std(SNR_indirect_baseline_raw')');
errorbar(p_list*1e3/50,mean(SNR_indirect_Lesica_raw,2),std(SNR_indirect_Lesica_raw')');
errorbar(p_list*1e3/50,mean(SNR_KL_full_raw,2),std(SNR_KL_full_raw')');
errorbar(p_list*1e3/50,mean(SNR_KL_woH_raw,2),std(SNR_KL_woH_raw')');
xlabel('Number of Background Noise Neurons (assuming each neuron is 50 Hz)');
legend('MISE','SNR indirect raw','SNR indirect baseline raw','SNR indirect Lesica raw', ...
    'SNR KL three componets','SNR KL two componets');
%title({'Mean Integrated Squared Error';'between True Stimulus and Normalized PSTH'});


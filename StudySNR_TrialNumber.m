clearvars;

addpath(genpath('D:/Github/GLM-RW'))

repnum = 3e3;
maxSig = 7e-2;
minSig = 2e-2;
signalType = 3; % 1 for no signal, 2 for square wave, 3 for gamma, 4 for Gaussian random
dt = 1e-3;
I_per = zeros(1,1e3);
if signalType == 2
    I_per(301:410) = 1;
end
if signalType == 3
    pd = makedist('InverseGaussian','mu',2e2,'lambda',6e2);
    I_per(201:1000) = pdf(pd,[1:800]);
end
if signalType == 4
    I_per = abs(normrnd(0.05,0.03,1,1e3));
end
if signalType == 5
    I_per = repmat([ones(1,50),zeros(1,50)],1,10);
end
I_per = I_per/max(I_per)*(maxSig-minSig) + minSig;
%plot(I_per);
I = repmat(I_per,repnum,1);
FRmulit = I;
ST = poissrnd(FRmulit);
FR = I_per;

meanFR = mean(FR);
SNR_direct = sum((FR-meanFR).^2*dt)/meanFR;
SNR_indirect = var(mean(ST)) / mean(var(ST)); 
SNR_indirect_baseline = var(mean(ST(:,201:800))) / mean(var(ST(:,[1:200,801:1000])));
SNR_indirect_Lesica = var(mean(ST)) / mean( var( (ST-mean(ST))' )' );
db_true = 10*log10(SNR_direct);
db_sq = 10*log10(SNR_indirect);

%% Plot Raster
N = 100;
T = 1e3;
y = ~(ST==0);

figure
subplot(2,1,1);
plotraster(y(1:N,1:T),1:T,'Simulated Result');
title('Spike train');
xlabel('ms');
ylabel('trial');

%% Tracking Signal
subplot(2,1,2);
hold on
tot_N = repnum;
fr = sum(ST);
ddt = 1;
if ddt ~= 1
    fr2 = reshape(fr,ddt,[]);
    fr2 = sum(fr2);
else
    fr2 = fr;
end
yyaxis left
semilogy( ((1:length(fr2))-1/2)*ddt ,fr2/tot_N/ddt );
ylim([0.01,0.08]);
%axis([0 5000 3 9]);
xlabel('t');
ylabel('Spike Probability');
yyaxis right
plot(I_per);
title('Spike Probability & Signal');
ylim([0.01,0.08]);
%axis([0 1000 -0.005 0.09]);
xlabel('t');
ylabel('Signal');

%% GLM with stimulus filter
y_glm = reshape(ST',1,[]);
y_glm(find(y_glm>1))=1;
I_AllInputs = repmat(FR,1,repnum);
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

%% without stimulus filter
fit_k = 0;

[k, h, dc, prs, kbasis, hbasis, stats] = ...
    fit_glm(zeros(1,length(I_AllInputs))',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
k = zeros(1,length(k))';
[pvalue, rate, h_output, k_output] = KStest(y_glm, h', I_AllInputs, k', dc, plotKS);
nlogl_0H = -sum( log(rate).*y_glm - rate );    % 0H : baseline, history included
Dev_0H = 2* (nlogl_0H - nlogl_satured);
%% without history filter
fit_k = 2;

[k, h, dc, prs, kbasis, hbasis, stats] = ...
    fit_glm(I_AllInputs',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
[pvalue, rate, h_output, k_output] = KStest(y_glm, h', I_AllInputs, k', dc, plotKS);
nlogl_0S = -sum( log(rate).*y_glm - rate );    % 0S : baseline, stimulus included
Dev_0S = 2* (nlogl_0S - nlogl_satured);

%% without both history filter and stimulus filter, only baseline
rate0 = ones(1,length(y_glm)) * sum(y_glm)/length(y_glm);
nlogl_0 = -sum( log(rate0).*y_glm - rate0 );    % 0 : baseline included
Dev_0 = 2* (nlogl_0 - nlogl_satured);

%% 
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

%% To plot SNR - trial number figure

recNeuNumList = ceil(logspace(0.5,2,10)); % for compare different varience-based
%recNeuNumList = [25,1e2,3e2,1e3,3e3];   % for compare KL and varience
repeat = 100;
db_raw = NaN*zeros(length(recNeuNumList),repeat);
db_rawL1 = NaN*zeros(length(recNeuNumList),repeat);
db_rawL2 = NaN*zeros(length(recNeuNumList),repeat);
db_baseline_raw = NaN*zeros(length(recNeuNumList),repeat);
db_Lesica_raw = NaN*zeros(length(recNeuNumList),repeat);

maxSig = 4e-2;  % 7e-2, 12e-2, 4e-2
minSig = 2e-2;
signalType = 3; % 1 for no signal, 2 for square wave, 3 for gamma, 4 for Gaussian random
I_per = zeros(1,1e3);
if signalType == 2
    I_per(301:410) = 1;
end
if signalType == 3
    pd = makedist('InverseGaussian','mu',2e2,'lambda',6e2);
    I_per(201:1000) = pdf(pd,[1:800]);
end
if signalType == 4
    I_per = abs(normrnd(0.05,0.03,1,1e3));
end
if signalType == 5
    I_per = repmat([ones(1,50),zeros(1,50)],1,10);
end
I_per = I_per/max(I_per)*(maxSig-minSig) + minSig;
%plot(I_per);

for ii = 1:length(recNeuNumList)
    for jj = 1:repeat
        dt = 1e-3;
        repnum = recNeuNumList(ii);

        I = repmat(I_per,repnum,1);
        FRmulit = I;
        ST = poissrnd(FRmulit);
        FR = I_per;

        meanFR = mean(FR);
        SNR_direct = sum((FR-meanFR).^2*dt)/meanFR;
        db_true = 10*log10(SNR_direct);
        SNR_indirect = var(mean(ST)) / mean(var(ST));
        SNR_indirect_baseline = var(mean(ST(:,201:800))) / mean(var(ST(:,[1:200,801:1000])));
        SNR_indirect_Lesica = var(mean(ST)) / mean( var( (ST-mean(ST))' )' );
        db_V = 10*log10(SNR_indirect);
        db_raw(ii,jj) = db_V;
        db_baseline_raw(ii,jj) = 10*log10(SNR_indirect_baseline);
        db_Lesica_raw(ii,jj) = 10*log10(SNR_indirect_Lesica);
%{
        % GLM with stimulus filter
        y_glm = reshape(ST',1,[]);
        y_glm(find(y_glm>1))=1;
        I_AllInputs = repmat(FR,1,repnum);
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
        
        db_rawL1(ii,jj) = db_glmS;
        db_rawL2(ii,jj) = db_glmS2;
        %}
    end
end

figure;
hold on
errorbar(recNeuNumList,mean(db_raw,2),std(db_raw')'/sqrt(repeat));
errorbar(recNeuNumList,mean(db_baseline_raw,2),std(db_baseline_raw')'/sqrt(repeat));
errorbar(recNeuNumList,mean(db_Lesica_raw,2),std(db_Lesica_raw')'/sqrt(repeat));
%errorbar(recNeuNumList,mean(db_rawL1,2),std(db_rawL1')');
%errorbar(recNeuNumList,mean(db_rawL2,2),std(db_rawL2')');
%plot([min(recNeuNumList),max(recNeuNumList)],[db_true,db_true]);
xlabel('Number of Trials');
ylabel('SNR /db');
title('SNR');
legend('Variance-based SNR','Variance-based SNR (baseline)',...
    'Variance-based SNR (Lesica)','Ground True SNR');
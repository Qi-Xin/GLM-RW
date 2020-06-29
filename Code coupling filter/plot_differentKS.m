clearvars;

addpath(genpath('D:/code'))

tau_E = 1e-3;           % 1ms
tau_I = tau_E;
tau_M = 20;
dt = 1;
p = 5e1;
q = 5e1;
bin = 5;   %ms
V_E = 0.023;
V_I = V_E;

plotFlag = 1;    % plot fit
inpNeuNum = 1e2;
tot_N = 1e4;
V_th = 1;
V_reset = 0;

signalType = 2; % 1 for no signal, 2 for square wave, 3 for gamma
I_per = 0*ones(1,1e3);
maxSig = 5e-2;
if signalType == 2
    I_per(101:400) = ones(1,300);
    I_per(601:800) = ones(1,200);
    I_per = I_per/max(I_per)*maxSig;
end
if signalType == 3
    pd = makedist('InverseGaussian','mu',2e2,'lambda',6e2);
    I_per(1:1000) = pdf(pd,[1:1000]);
    I_per = I_per/max(I_per)*maxSig;
end


%recNeuNumList = [1,5,20,50,99];
recNeuNumList = [1,10,100];
sessionList = [10,30,60];
repeat = 100;
IMSE_raw = NaN*zeros(length(recNeuNumList),repeat);
nlogl_raw = NaN*zeros(length(recNeuNumList),repeat);
IMSE_fix1 = NaN*zeros(1,repeat);
nlogl_fix1 = NaN*zeros(1,repeat);

for ii = 1:length(recNeuNumList)
    for jj = 1:length(sessionList)
        for kk = 1:repeat
            recNeuNum = recNeuNumList(ii);
            session = sessionList(jj);

            tot_t = session*1e3;
            repNum = session;
            I = repmat(I_per,1,repNum);
            I = I(1:tot_t);
            I_noInp = zeros(1,tot_t);


            I_AllInputsInd = random('poisson',repmat(I,inpNeuNum,1),inpNeuNum,tot_t);
            I_AllInputsInd(find(I_AllInputsInd>=2)) = 1;
            I_AllInputs = sum(I_AllInputsInd);
            I_AllInputs = I_AllInputs/inpNeuNum;

            randomChoose = randperm(inpNeuNum);
            randomChoose = randomChoose(1:recNeuNum);
            I_record = sum(I_AllInputsInd(randomChoose,:),1);


            % Simulation
            [ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I_AllInputs,tot_t,dt);
            y = full(y_sparse);


            % Fit Pillow GLM (Recorded Neuron)
            T = tot_t;
            %T = 1e6;
            y_glm = y(1,1:T);

            nkt = 150; % number of ms in stim filter
            kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
            kbasprs.ncos = 7; % number of raised-cosine vectors to use
            kbasprs.kpeaks = [1 round(nkt/1.5)];  % position of first and last bump (relative to identity bumps)
            kbasprs.b = 5; % how nonlinear to make spacings (larger -> more linear)
            %%% basis functions for post-spike kernel
            ihbasprs.ncols = 10;  % number of basis vectors for post-spike kernel
            hPeaksMax = 50;
            ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
            ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
            ihbasprs.absref = 0; % absolute refractory period, in ms

            fit_k = 1;
            if max(I) == 0
                fit_k = 0;
            end

            plotKS = 0;

            [kn, hn, dcn, prsn, kbasis, hbasis, stats] = fit_glm(I_record',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,0);
            [pvaluen, raten, h_outputn, k_outputn] = KStest(y_glm, hn', I_record, kn', dcn, plotKS);
            nlogln = -sum(log( raten.*y_glm + (1-raten).*(1-y_glm) ));

            p_raw(ii,jj,kk) = pvaluen;
        end
    end
end

figure;
hold on
for ii = 1:length(recNeuNumList)
    errorbar(sessionList+ii/1e1,mean(p_raw(ii,:,:),3),std(p_raw(ii,:,:),0,3)./sqrt(repeat));
end
xlabel('Number of Neurons recorded');
ylabel('p value');
title('p value');


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
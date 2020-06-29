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

inpNeuNum = 1e2;
recNeuNum = 1;
tot_t = 3e6;
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

repNum = ceil(tot_t/length(I_per));
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


%% Simulation
[ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I_AllInputs,tot_t,dt);
y = full(y_sparse);

%% Plot Raster
T = 1e3;
N = 100;

figure
subplot(2,1,1);
plotraster(reshape(y(1:N*T),[],N)',1:T,'Simulated Result');
title('Spike train');
xlabel('ms');
ylabel('trial');


%% Tracking Signal
ddt = bin;
subplot(2,1,2);
hold on
tot_N = (tot_t/length(I_per));
fr = sum(reshape(y,[],tot_N)');
if ddt ~= 1
    fr2 = reshape(fr,ddt,[]);
    fr2 = sum(fr2);
else
    fr2 = fr;
end
yyaxis left
semilogy( ((1:length(fr2))-1/2)*ddt ,fr2/tot_N/ddt );
%axis([0 5000 3 9]);
xlabel('t');
ylabel('Spike Probability');
yyaxis right
plot(I_per);
title('Spike Probability & Signal');
%axis([0 1000 -0.005 0.09]);
xlabel('t');
ylabel('Signal');
%% Fit Pillow GLM (Recorded Neuron)
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
plotFlag = 1;    % plot fit
plotKS = 0;

[kn, hn, dcn, prsn, kbasis, hbasis, stats] = fit_glm(I_record',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
[pvaluen, raten, h_outputn, k_outputn] = KStest(y_glm, hn', I_record, kn', dcn, plotKS);
nlogln = -sum(log( raten.*y_glm + (1-raten).*(1-y_glm) ));

%% Fit Pillow GLM (Population) 
T = tot_t;
%T = 1e6;
y_glm = y(1,1:T);

nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 5; % number of raised-cosine vectors to use
kbasprs.kpeaks = [1 round(nkt/1.5)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 5; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 10;  % number of basis vectors for post-spike kernel
hPeaksMax = 50;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

[kp, hp, dcp, prsp, kbasis, hbasis, stats] = fit_glm(I_AllInputs',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
[pvaluep, ratep, h_outputp, k_outputp] = KStest(y_glm, hp', I_AllInputs, kp', dcp, plotKS);
nloglp = -sum(log( ratep.*y_glm + (1-ratep).*(1-y_glm) ));
%% Plot Output
T = 1e3;
figure
subplot(3,1,1)
hold on
plot(1:T,h_outputn(1:T),'b');
plot(1:T,h_outputp(1:T),'r');
title('Post-spike Filter Output');
legend('Coupling Filter','Stimulus Filter');

subplot(3,1,2)
hold on
plot(1:T,k_outputn(1:T),'b');
plot(1:T,k_outputp(1:T),'r');
N = 1;
yyaxis right
plotraster(reshape(I_record(1:N*T),[],N)',1:T,'');
ylim([-10,1]);
title('Stilmulus / Coulping Filter Output');
legend('Coupling Filter','Stimulus Filter');

subplot(3,1,3)
yyaxis left
hold on
plot(1:T,raten(1:T),'b');
plot(1:T,ratep(1:T),'r-');
title('Predicted Firing Rate');
legend('Coupling Filter','Stimulus Filter');

N = 1;
yyaxis right
plotraster(reshape(y(1:N*T),[],N)',1:T,'');
ylim([-10,1]);

%% Fit Hazard
%{

nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 5; % number of raised-cosine vectors to use
kbasprs.kpeaks = [1 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 15;  % number of basis vectors for post-spike kernel
hPeaksMax = 100;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

[k, h, dc, prs, kbasis, hbasis , stats] = fit_hazard(I_record',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
%}
%% Fit FNL GLM 
%{
NumSP = 3;
T = tot_t;
%T = 1e6;
y_glm = y(1,1:T);

nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 5; % number of raised-cosine vectors to use
kbasprs.kpeaks = [1 round(nkt/1.5)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 5;  % number of basis vectors for post-spike kernel
hPeaksMax = 70;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 1*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

fit_k = 1;
if max(I) == 0
    fit_k = 0;
end
plotFlag = 1;    % plot fit

[k, h, dc, prs, kbasis, hbasis, stats] = fit_FNL(I',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag,NumSP,ISI);
%}
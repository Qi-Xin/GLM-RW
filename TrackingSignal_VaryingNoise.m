clearvars;

tau_E = 1e-3;           % 1ms
tau_I = tau_E;
tau_M = 20;
dt = 1;
p = 5e1;
q = 5e1;
%{
p = 0.3e1;
q = 0e1;

p = 5e1;
q = 5e1;

p = 5e-2;
q = 5e-2;

p = 1e4;
q = 1e4;
%}

tot_t = 1e7;
bin = 1;   %ms
ddt = bin;
V_E = 0.02;
V_I = 0.02;
adjStepOrNot = 0;
adjValue = 50;
V_th = 1;
V_reset = 0;

maxSig = 0.9e-1;
signalType = 3; % 1 for no signal, 2 for square wave, 3 for gamma
I_per = zeros(1,1e3);
if signalType == 2
    %I_per(301:410) = 1;
    I_per(1:1000) = 1;
    %{
    I_per(101:110) = 0.1;
    I_per(201:210) = 0.2;
    I_per(301:310) = 0.3;
    I_per(401:410) = 0.4;
    I_per(501:510) = 0.5;
    I_per(601:610) = 0.6;
    I_per(701:710) = 0.7;
    I_per(801:810) = 0.8;
    I_per(901:910) = 0.9;
    %}
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
%I = zeros(T,1);

%{
Inverse Gaussian: 
tau_E = 10*1e-5;           % 0.1ms
tau_I = tau_E;          % 0.1ms
tau_M = 10*1e5;          % 0.1ms
dt = 10*0.1;            % 0.1ms
p = 1.2e-2;
q = 1e-2;

Clock-like: 
tau_E = 10*1e-5;           % 0.1ms
tau_I = tau_E;          % 0.1ms
tau_M = 10*500;          % 0.1ms
dt = 10*0.1;            % 0.1ms
p = 3e-1;
q = 1e-1;

Poission-like: 
tau_E = 10*1e-5;           % 0.1ms
tau_I = tau_E;          % 0.1ms
tau_M = 10*25;          % 0.1ms
dt = 10*0.1;            % 0.1ms
p = 1.2e-2;
q = 1e-2;

Burst: 
tau_E = 10*2e2;           % 0.1ms
tau_I = tau_E;          % 0.1ms
tau_M = 10*25;          % 0.1ms
dt = 10*0.1;            % 0.1ms
p = 0.5e-4;
q = 0.2e-4;

tau_E = 10*1e-3;           % 0.1ms   for illustrating tau_input: 1e-3,30,200
tau_I = tau_E;          % 0.1ms
tau_M = 10*25;          % 0.1ms
dt = 10*0.1;            % 0.1ms
p = 2e-2;
q = 1e-2;
%}


%% Adjust Steps
if adjStepOrNot == 1
    x_up = 100;
    x_down = 1e-5;
    error = 1e-5;
    res_down = GetMeanISI_J(tau_E,tau_I,tau_M,x_down*V_E,x_down*V_I,p,q,V_th,V_reset,I_noInp,tot_t,dt,adjValue);
    res_up = GetMeanISI_J(tau_E,tau_I,tau_M,x_up*V_E,x_up*V_I,p,q,V_th,V_reset,I_noInp,tot_t,dt,adjValue);
    while(res_down * res_up < 0)
        x = 0.5*(x_up + x_down);
        res = GetMeanISI_J(tau_E,tau_I,tau_M,x*V_E,x*V_I,p,q,V_th,V_reset,I_noInp,tot_t,dt,adjValue);
        if( res*res_down < 0 )
            x_up = x;
        else
            x_down = x;
        end
        if( abs(x_up-x_down) < error )
            break;
        end
    end
    result_x = 0.5*(x_up + x_down)
    V_E = result_x*V_E;
    V_I = result_x*V_I;
end
%% Simulation
[ISI,spike_timing,y_sparse,V,inputE,inputI,p_varying,q_varying] = ...
GetISI_VaryingNoise(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt,[16,16],[-1,0.1]);
%{
[100,100],[0,0]

[100,100],[0.3,0.3]
[100,100],[0,0.01]

[5,5],[-1,0.1]
[5,5],[-1,0.01]
%}
y = full(y_sparse);

figure
subplot(2,1,1);
plotraster(reshape(y(1:N*T),[],N)',1:T,'Simulated Result');
title('Spike train');
xlabel('ms');
ylabel('trial');

% Tracking Signal
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
hold on
plot(I_per/max(I_per));
plot(q_varying(1:1e3)/mean(q_varying));
plot(p_varying(1:1e3)/mean(p_varying));
legend("PSTH","Signal","E noise","I noise")
%axis([0 1000 -0.005 0.09]);
xlabel('t');
ylabel('Signal & Noise Amplitute');
title('Spike Probability ');
yyaxis left
%% MISE
target = fr2/tot_N/ddt;
groundTrue = I_per;
%baseline = mean(target(1:200));
%fun = @(k) sum(((target-baseline)*k - groundTrue).^2)
%x0 = 1;
fun = @(k) sum(((target-k(2))*k(1) - groundTrue).^2)
x0 = [1,0.2];
[x,fval] = fminunc(fun,x0)

%% Fit GLM with no input and no input filter

T = tot_t;
%T = 1e6;
y_glm = y(1,1:T);

nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 5; % number of raised-cosine vectors to use
kbasprs.kpeaks = [1 round(nkt/1.5)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 10;  % number of basis vectors for post-spike kernel
hPeaksMax = 100;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

fit_k = 1;
if max(I) == 0
    fit_k = 0;
end
plotFlag = 1;    % plot fit

[kn, hn, dcn, prsn, kbasis, hbasis, stats] = fit_glm(I',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
[pvaluen, raten, h_outputn, k_outputn] = KStest(y_glm, hbasis(:,3)', I, kbasis(:,3)', dcn, 0);
figure
hold on
T = 1e3;
plot(1:T,h_outputn(1:T),'b');
plot(1:T,k_outputn(1:T),'r');
title('A Filter Basis Output');
legend('History','Stimulus');
%% Simulation for baseline hazard
%{
tot_t = 1e6;
I_noInp = zeros(1,tot_t);
[ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I_noInp,tot_t,dt);
y = full(y_sparse);

% Fit Hazard
T = tot_t;
%T = 1e6;
y_glm = y(1,1:T);
fit_k = 0;
plotFlag = 1;

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

figure
[k, h, dc, prs, kbasis, hbasis , stats] = fit_hazard(I_noInp',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
%}
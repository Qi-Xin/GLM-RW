clearvars;

tau_E = 1e-3;           % 1ms
tau_I = tau_E;
tau_M = 20;
dt = 1;
p = 5e1;
q = 5e1;

bin = 5;   %ms
%V_E = 1*(1-exp(-dt/tau_E));
V_E = 0.023;
%V_E = 0.1;
V_I = V_E;
adjStepOrNot = 0;
adjValue = 50;        % 0.1ms
tot_t = 1e6;
tot_N = 1e4;
V_th = 1;
V_reset = 0;

signalType = 2; % 1 for no signal, 2 for square wave, 3 for gamma
I_per = 0*ones(1,1e3);
maxSig = 5e-2;
if signalType == 2
    I_per(251:750) = ones(1,500);
    I_per = I_per/max(I_per)*maxSig;
end
if signalType == 3
    pd = makedist('InverseGaussian','mu',2e2,'lambda',6e2);
    I_per(1:1000) = pdf(pd,[1:1000]);
    I_per = I_per/max(I_per)*maxSig;
end
%plot(I_per);


repnum = ceil(tot_t/length(I_per));
I = repmat(I_per,1,repnum);
I = I(1:tot_t);
I_noInp = zeros(1,tot_t);

I_record = random('poisson',I,1,tot_t);
I_record(find(I_record>=2)) = 1;
I_record = maxSig*I_record;
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
[ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt);
y = full(y_sparse);
%% 
ddt = bin;
max1 = ceil(max(ISI));
dist = zeros(1,ceil(max1/ddt));
total_trial = length(ISI);
for j = 1:total_trial
    dist(ceil(ISI(j)/ddt)) = dist(ceil(ISI(j)/ddt))+1;
end

dead = 0;
svv(1) = total_trial;
for i = 1:length(dist)
    dead = dead + dist(i);
    svv(i+1) = total_trial-dead;
end
svv = svv(1:length(svv)-1);
[rate,err] = binofit(dist,svv);
xx1 = ((1:length(dist))-1/2);
xdata = (xx1)*ddt;
ydata = rate;
last = max(find(diff(err')<=0.1));
xdata = xdata(1:last);
ydata = ydata(1:last);

figure
subplot(3,1,1)
yyaxis left
h = histogram(ISI,0:ddt:max1,'Normalization','pdf');
xlabel('t/ms');
ylabel('ISI distribution');
xlim([0 200]);
%axis([0 1000 0 0.0075]);
yyaxis right
errorbar(xdata,rate(1:last)/ddt,err(1:last,1)'/ddt-rate(1:last)/ddt,-err(1:last,2)'/ddt+rate(1:last)/ddt);
ylabel('Hazard function');
xlim([0 200])
%axis([0 1000 0 0.02]);
title(['tau_E = ',num2str(tau_E), ...
    ';tau_I = ',num2str(tau_I), ...
    ';tau_M = ',num2str(tau_M), ...
    ';V_E = ',num2str(V_E), ...
    ';V_I = ',num2str(V_I), ...
    ';p = ',num2str(p), ...
    ';q = ',num2str(q)]);

%{
title('Hazard function');
legend('r2');
xlabel('t/s');
ylabel('rate');
axis([0 200 0 1]);
%}

%% Plot Raster
T = 3e3;
N = 10;
%{
figure
plotraster(reshape(y(1:10*T),[],10)',1:T,'Simulated Result');
%}
subplot(6,1,3);
plot((1:T),inputE(1:T),(1:T),inputI(1:T),(1:T),I(1:T));
xlabel('t/s');
ylabel('Amplitude');
legend('Excitatory Poisson Input','Inhibitory Poisson Input');
title('Poisson Noise Input');
subplot(6,1,4);
hold on
plot((1:T),I(1:T));
plot((1:T),I_record(1:T));
xlabel('t/s');
ylabel('Amplitude');
legend('Population Input','Recorded Neuron');
title('Signal Input');
ylim([0 1e-1]);
subplot(6,1,5);
plot(1e-4*(1:T),V(1:T));
xlabel('t/s');
ylabel('Voltage');title('Voltage');
subplot(6,1,6);
plotraster(reshape(y_sparse(1:1*T),[],1)',1:T,'Spike train');
title('Spike train');
%% Plot Raster
T = 1e3;
N = 100;

figure
subplot(2,1,1);
plotraster(reshape(y(1:N*T),[],N)',1:T,'Simulated Result');
title('Spike train');
xlabel('ms');
ylabel('trial');

%{
plotraster(reshape(y(1:10*T),[],10)',1:T,'Simulated Result');

subplot(6,1,3);
plot((1:T),inputE(1:T),(1:T),inputI(1:T),(1:T),I(1:T));
xlabel('t/s');
ylabel('Amplitude');
legend('Excitatory Poisson Input','Inhibitory Poisson Input');
title('Poisson Noise Input');
subplot(6,1,4);
plot((1:T),I(1:T));
xlabel('t/s');
ylabel('Amplitude');
legend('Signal Input');
title('Signal Input');
ylim([0 1e-1]);
subplot(6,1,5);
plot(1e-4*(1:T),V(1:T));
xlabel('t/s');
ylabel('Voltage');
title('Voltage');
subplot(6,1,6);
plotraster(reshape(y_sparse(1:1*T),[],1)',1:T,'Spike train');
title('Spike train');
%}
%% Tracking Signal
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

nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 5; % number of raised-cosine vectors to use
kbasprs.kpeaks = [1 round(nkt/1.5)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 5;  % number of basis vectors for post-spike kernel
hPeaksMax = 100;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

fit_k = 1;
if max(I) == 0
    fit_k = 0;
end
plotFlag = 1;    % plot fit

[kn, h, dc, prs, kbasis, hbasis, stats] = fit_glm(I_record',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
%% Fit Pillow GLM (Population) 
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
hPeaksMax = 100;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

fit_k = 1;
if max(I) == 0
    fit_k = 0;
end
plotFlag = 1;    % plot fit

[kp, h, dc, prs, kbasis, hbasis, stats] = fit_glm(I',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
%%
T=1e4;
outn=sameconv(I_record',flipud(kn));
outp=sameconv(I',flipud(kp));
figure
hold on
plot(1:T,outn(1:T));
plot(1:T,outp(1:T))
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
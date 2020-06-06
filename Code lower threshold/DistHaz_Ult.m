clearvars;

tau_E = 10*50;           % 0.1ms
tau_I = tau_E;          % 0.1ms
tau_M = 10*50;          % 0.1ms
dt = 10*0.1;            % 0.1ms
p = 6e-2;
q = 3e-2;
th = 1;
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

bin = 10;   %ms
V_E = 1*(1-exp(-dt/tau_E));
V_I = 1*(1-exp(-dt/tau_I));
adjStepOrNot = 0;
adjValue = 10*50;        % 0.1ms
tot_t = 1e6;
tot_N = 1e4;
V_th = 1;
V_reset = 0;

%% Adjust Steps
if adjStepOrNot == 1
    x_up = 100;
    x_down = 1e-5;
    error = 1e-5;
    res_down = GetMeanISI_J(tau_E,tau_I,tau_M,x_down*V_E,x_down*V_I,p,q,V_th,V_reset,tot_N,tot_t,dt,adjValue);
    res_up = GetMeanISI_J(tau_E,tau_I,tau_M,x_up*V_E,x_up*V_I,p,q,V_th,V_reset,tot_N,tot_t,dt,adjValue);
    while(res_down * res_up < 0)
        x = 0.5*(x_up + x_down);
        res = GetMeanISI_J(tau_E,tau_I,tau_M,x*V_E,x*V_I,p,q,V_th,V_reset,tot_N,tot_t,dt,adjValue);
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
[ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI(th,tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,tot_N,tot_t,dt);
ISI = 0.1*ISI;
spike_timing = 0.1*spike_timing;
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
%{
figure
subplot(2,1,1)
yyaxis left
h = histogram(ISI,0:ddt:max1,'Normalization','pdf');
xlabel('t/ms');
ylabel('ISI distribution');
%axis([0 1000 0 0.0075]);
yyaxis right
errorbar(xdata,rate(1:last)/ddt,err(1:last,1)'/ddt-rate(1:last)/ddt,-err(1:last,2)'/ddt+rate(1:last)/ddt);
ylabel('Hazard function');
xlim([0 1000])
%axis([0 1000 0 0.02]);
title(['tau_E = ',num2str(tau_E), ...
    ';tau_I = ',num2str(tau_I), ...
    ';tau_M = ',num2str(tau_M), ...
    ';V_E = ',num2str(V_E), ...
    ';V_I = ',num2str(V_I), ...
    ';p = ',num2str(p), ...
    ';q = ',num2str(q)]);
%}
%{
title('Hazard function');
legend('r2');
xlabel('t/s');
ylabel('rate');
axis([0 200 0 1]);
%}

%% Plot Raster
T = (10*1e4);
N = 10;
%{
figure
plotraster(reshape(y(1:10*T),[],10)',1:T,'Simulated Result');
%}
%{
subplot(6,1,4);
plot(1e-4*(1:T/1),inputE(1:T/1)-inputI(1:T/1));
xlabel('t/s');
ylabel('Amplitude');
title('Equivalent Poisson Input');
subplot(6,1,5);
plot(1e-4*(1:T),V(1:T));
xlabel('t/s');
ylabel('Normalized Voltage');
title('Voltage');
subplot(6,1,6);
plotraster(reshape(y_sparse(1:1*T),[],1)',1:T,'Spike train');
%}

%% Fit GLM with no input and no input filter
%{
T = tot_t;
%T = 1e6;
I = zeros(T,1);
y_glm = y(1,1:T);

nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 3; % number of raised-cosine vectors to use
kbasprs.kpeaks = [5 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 10;  % number of basis vectors for post-spike kernel
hPeaksMax = 3*max([mean(ISI) tau_E/10 tau_I/10 tau_M/10]);
ihbasprs.hpeaks = [5 10*hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 10*0.5*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

fit_k = 1;
if max(I) == 0
    fit_k = 0;
end
plotFlag = 1;    % plot fit

[k, h, dc, prs, kbasis, hbasis] = fit_glm(I,y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
%}
%% Fit Hazard
%{
nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 3; % number of raised-cosine vectors to use
kbasprs.kpeaks = [5 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 15;  % number of basis vectors for post-spike kernel
hPeaksMax = 10*max([mean(ISI) tau_E/10 tau_I/10 tau_M/10]);
ihbasprs.hpeaks = [1 10*hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 10*0.5*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

[k, h, dc, prs, kbasis, hbasis] = fit_hazard(I,y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
%}
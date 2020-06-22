clearvars;

tot_t = 2e6;
%tot_N = 1e4;
V_th = 1;
V_reset = 0;
dt = 1;

signalType = 1; % 1 for no signal, 2 for square wave, 3 for gamma
I_per = zeros(1,1e3);
maxSig = 5e-2;
if signalType == 2
    I_per(251:750) = ones(1,500);
    I_per = I_per/max(I_per)*maxSig;
end
I_noInp = zeros(1,tot_t);
if signalType == 3
    pd = makedist('InverseGaussian','mu',2e2,'lambda',6e2);
    I_per(201:1000) = pdf(pd,[1:800]);
    I_per = I_per/max(I_per)*maxSig;
end
repnum = ceil(tot_t/length(I_per));
I = repmat(I_per,1,repnum);
I = I(1:tot_t);
I_noInp = zeros(1,tot_t);

e_range = linspace(1e-1,1.5e1,10);
mask_range = [0,1];
V_E = 1.5e-1;
V_I = 1e-1;
tau_M = 20;
tau_E = 10;
tau_I = 10;

adjStepOrNot = 0;
adjValue = 50;
le = length(e_range);
lm = length(mask_range);
hazard_rcd = cell(le,lm);
history_rcd = cell(le,lm);
stimulus_rcd = cell(le,lm);
hazardLU_rcd = cell(le,lm);
historyLU_rcd = cell(le,lm);
stimulusLU_rcd = cell(le,lm);
bias_rcd = cell(le,lm);
biasLU_rcd = cell(le,lm);
y_rcd = cell(le,lm);
y_sparse_rcd = cell(le,lm);
ISI_rcd = cell(le,lm);
spike_timing_rcd = cell(le,lm);
mean_rcd = NaN*zeros(le,lm);
cv_rcd = NaN*zeros(le,lm);

rateBias = -3;
w = -20*exp(-(1:50)/3);
w = [-20*ones(1,20),w];
simulationtimes = 1;
T = 1e5;
nHistBins = numel(w);
for i = 1:(le)
    for j = 1:(lm)
        [i,j]
        
        p = e_range(i);
        q = e_range(i) + mask_range(j);
        
        [~,~,~,~,Ie,Ii] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,T+nHistBins,dt);
        
        y = zeros(simulationtimes , T + nHistBins);
        for NSimTrial=1:simulationtimes
            for t = (nHistBins+1):(T+nHistBins)
                yy = poissrnd(exp(w * (y(NSimTrial,t - (1:nHistBins)))' + rateBias + Ie(t)-Ii(t)));
                if yy ~= 0
                    y(NSimTrial,t) = 1;
                end
            end

        end
        y = y(:,nHistBins+1:end);
        fr(i,j) = sum(y)/T*1e3;
    end
end
        
        %{
        if signalType ~= 1
            % GLM
            ISI = ISI_rcd{i,j};
            y = y_rcd{i,j};
            T = tot_t;
            %I = zeros(T,1);
            y_glm = y(1,1:T);
            
            nkt = 100; % number of ms in stim filter
            kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
            kbasprs.ncos = 5; % number of raised-cosine vectors to use
            kbasprs.kpeaks = [1 round(nkt/1.5)];  % position of first and last bump (relative to identity bumps)
            kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
            %%% basis functions for post-spike kernel
            ihbasprs.ncols = 7;  % number of basis vectors for post-spike kernel
            hPeaksMax = 100;
            ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
            ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
            ihbasprs.absref = 0; % absolute refractory period, in ms
            
            fit_k = 1;
            if max(I) == 0
                fit_k = 0;
            end
            plotFlag = 0;    % plot fit
            
            [k, h, dc, prs, kbasis, hbasis, stats, kLU, hLU, dcLU] = fit_glm(I',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
            stimulus_rcd{i,j} = k;
            stimulusLU_rcd{i,j} = kLU;
            history_rcd{i,j} = h;
            historyLU_rcd{i,j} = hLU;
            bias_rcd{i,j} = dc;
            biasLU_rcd{i,j} = dcLU;
            
        else
            % hazard
            ihbasprs.ncols = 7;  % number of basis vectors for post-spike kernel
            hPeaksMax = 50;
            ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
            ihbasprs.b = 0.1*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
            ihbasprs.absref = 0; % absolute refractory period, in ms
            
            [k, h, dc, prs, kbasis, hbasis, stats, kLU, hLU] = fit_hazard(I',y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag); 
            hazard_rcd{i,j} = h;
            hazardLU_rcd{i,j} = hLU;
        end
        %}

xdata = e_range;
ydata = fr(:,1)';
r = @(k,xdata) k(1)./(k(2)+1./exp(k(3).*log(xdata)));
k0 = [100,5,5];
kfit1 = lsqcurvefit(r,k0,xdata,ydata);
XX = 0:1e-1:15;
YY1 = r(kfit1,XX);
%  5.2930    0.0497    1.2600

xdata = e_range;
ydata = fr(:,2)';
r = @(k,xdata) k(1)./(k(2)+1./exp(k(3).*log(xdata)));
k0 = [100,5,5];
kfit2 = lsqcurvefit(r,k0,xdata,ydata);
XX = 0:1e-1:15;
YY2 = r(kfit2,XX);
% 1.1790    0.0113    1.5953

figure;
hold on
plot(e_range,fr(:,1)',e_range,fr(:,2)');
plot(XX,YY1);
plot(XX,YY2);
legend('mask off','mask on');
xlabel('Input Firing Rate');
ylabel('Response Firing Rate');


%%
%{
%% Plot CV in colormap
cv_rcd_pcolor = NaN*zeros(size(cv_rcd)+1);
cv_rcd_pcolor(1:end-1,1:end-1) = cv_rcd;
pcolor(cv_rcd_pcolor);
colorbar;
%% Plot CV in lines
figure
%hold on
%semilogx(tau_E_range(2:end),cv_rcd(2:end,1));
%semilogx(tau_E_range(2:end),cv_rcd(2:end,2));
%semilogx(tau_E_range(2:end),cv_rcd(2:end,3));
semilogx(tau_E_range(2:end),cv_rcd(2:end,4));
legend('\tau_m=10','\tau_m=25','\tau_m=50','\tau_m=100');
xlabel('\tau_{input}');
ylabel('CV');
%% 
ISI = ISI_rcd{1,end};

bin = 2;
ddt = bin;
max1 = ceil(max(ISI));
dist = zeros(1,ceil(max1/ddt));
total_trial = length(ISI);
for j = 1:total_trial
    dist(ceil(ISI(j)/ddt)) = dist(ceil(ISI(j)/ddt))+1;
end

clear svv
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
yyaxis left
h = histogram(ISI,0:ddt:max1,'Normalization','pdf');
xlabel('t/ms');
ylabel('ISI distribution');
%axis([0 1000 0 0.0075]);
yyaxis right
errorbar(xdata,rate(1:last)/ddt,err(1:last,1)'/ddt-rate(1:last)/ddt,-err(1:last,2)'/ddt+rate(1:last)/ddt);
ylabel('Hazard function');
%axis([0 1000 0 0.015]);
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

% Plot Raster
%{
y = y_rcd{1,end};
T = (10*1e3);
N = 10;
figure
plotraster(reshape(y(1:10*T),[],10)',1:T,'Simulated Result');
%}
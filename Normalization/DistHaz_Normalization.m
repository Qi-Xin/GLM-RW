clearvars;

tot_t = 1e6;
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

e_range = logspace(-1,2,10);
mask_range = [0,10];
V_E = 0.1;
V_I = 0.1;
tau_M = 20;
tau_E = 1e-5;
tau_I = 1e-5;

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

for i = 1:(le)
    for j = 1:(lm)
        [i,j]
        
        p = e_range(i);
        q = e_range(i) + mask_range(j);
        % Simulation
        [ISI_rcd{i,j},spike_timing_rcd{i,j},y_sparse_rcd{i,j}] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt);
        mean_rcd(i,j) = mean(ISI_rcd{i,j});
        cv_rcd(i,j) = std(ISI_rcd{i,j})/mean_rcd(i,j);
        y_rcd{i,j} = full(y_sparse_rcd{i,j});
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
    end
end
fr = 1e3./mean_rcd;

figure;
semilogx(e_range,fr(:,1)',e_range,fr(:,2)');
legend('mask off','mask on');
xlabel('Input Firing Rate');
ylabel('Response Firing Rate');
%% Plot filter & hazard function
figure
nn = 0;
for i = 1:length(tau_E_range)
    for j = 1:length(tau_M_range)
        nn = nn + 1;
        subplot(length(tau_E_range),length(tau_M_range),nn);hold on
        h = history_rcd{i,j};
        hL = historyLU_rcd{i,j}(:,1);
        hU = historyLU_rcd{i,j}(:,2);
        plot((1:length(h)),h);
        fill([(1:length(h)) fliplr(1:length(h))],[hL' fliplr(hU')],'b','facealpha',0.2,'edgealpha',0);
        xlim([1 length(h) ]);
        title(['\tau_{input} = ',num2str(tau_E_range(i)),';\tau_M = ',num2str(tau_M_range(j))]);
        grid on
    end
end
figure
nn = 0;
for i = 1:length(tau_E_range)
    for j = 1:length(tau_M_range)
        nn = nn + 1;
        subplot(length(tau_E_range),length(tau_M_range),nn);hold on
        if fit_k == 1
            k = stimulus_rcd{i,j};
            kL = stimulusLU_rcd{i,j}(:,1);
            kU = stimulusLU_rcd{i,j}(:,2);
            plot((1:length(k)),k);
            fill([(1:length(k)) fliplr(1:length(k))],[kL' fliplr(kU')],'b','facealpha',0.2,'edgealpha',0);
            xlim([0 length(k)]);
            set(gca,'xtick',0:round(length(k)/2):length(k),'xticklabel',round(fliplr(-dt*(0:round(length(k)/2):length(k)))));
            grid on
            title(['\tau_{input} = ',num2str(tau_E_range(i)),';\tau_M = ',num2str(tau_M_range(j))]);
        else
            h = hazard_rcd{i,j};
            hL = hazardLU_rcd{i,j}(:,1);
            hU = hazardLU_rcd{i,j}(:,2);
            plot((1:length(h)),h);
            fill([(1:length(h)) fliplr(1:length(h))],[hL' fliplr(hU')],'b','facealpha',0.2,'edgealpha',0);
            xlim([1 length(h) ]);
            title(['\tau_{input} = ',num2str(tau_E_range(i)),';\tau_M = ',num2str(tau_M_range(j))]);
            grid on
        end

    end
end
figure
T = 1e3;
N = 10;
nn = 0;
for i = 1:length(tau_E_range)
    for j = 1:length(tau_M_range)
        nn = nn + 1;
        subplot(length(tau_E_range),length(tau_M_range),nn);
        plotraster(reshape(y_sparse_rcd{i,j}(T+1:N*T+T),[],N)',1:T,'');
        title(['\tau_{input} = ',num2str(tau_E_range(i)),';\tau_M = ',num2str(tau_M_range(j))]);
    end
end
%% Features
[peakVal,peakLoc,dipVal,dipLoc] = get_features(history_rcd,historyLU_rcd,0,tau_E_range,tau_M_range);
[peakVal,peakLoc,dipVal,dipLoc] = get_features(stimulus_rcd,stimulusLU_rcd,1,tau_E_range,tau_M_range);

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
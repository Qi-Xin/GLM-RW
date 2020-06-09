clearvars;

%tau_E_range = [1,5,10,20,50,100,200,500];
%tau_M_range = [10,25,50,100];
tau_E_range = [1,10,50,150,500];
tau_M_range = [10,25,50];
adjStepOrNot = 1;
adjValue = 10*50;
le = length(tau_E_range);
lm = length(tau_M_range);
hazard_rcd = cell(le,lm);
history_rcd = cell(le,lm);
bias_rcd = cell(le,lm);
y_rcd = cell(le,lm);
y_sparese_rcd = cell(le,lm);
ISI_rcd = cell(le,lm);
spike_timing_rcd = cell(le,lm);
mean_rcd = NaN*zeros(le,lm);
cv_rcd = NaN*zeros(le,lm);
for i = 1:length(tau_E_range)
    for j = 1:length(tau_M_range)
        display([i,j]);
        tau_E = 10*tau_E_range(i);
        tau_I = tau_E;
        tau_M = 10*tau_M_range(j);
        dt = 10*0.1;
        p = 6e-2;
        q = 3e-2;

        bin = 10;   %ms
        V_E = 1*(1-exp(-dt/tau_E));
        V_I = 1*(1-exp(-dt/tau_I));
        tot_t = 1e6;
        Nn = 1e3;
        tot_N = 1e4;
        V_th = 1;
        V_reset = 0;
        
        % Adjust step
        if adjStepOrNot == 1
            x_up = 100;
            x_down = 1e-5;
            error = 1e-5;
            res_down = GetMeanISI_J(tau_E,tau_I,tau_M,x_down*V_E,x_down*V_I,p,q,V_th,V_reset,tot_N,tot_t,dt,Nn,adjValue);
            res_up = GetMeanISI_J(tau_E,tau_I,tau_M,x_up*V_E,x_up*V_I,p,q,V_th,V_reset,tot_N,tot_t,dt,Nn,adjValue);
            while(res_down * res_up < 0)
                x = 0.5*(x_up + x_down);
                res = GetMeanISI_J(tau_E,tau_I,tau_M,x*V_E,x*V_I,p,q,V_th,V_reset,tot_N,tot_t,dt,Nn,adjValue);
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

        % Simulation
        [ISI_rcd{i,j},spike_timing_rcd{i,j},y_sparse_rcd{i,j}] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,tot_N,tot_t,dt,Nn);
        ISI_rcd{i,j} = 0.1*ISI_rcd{i,j};
        spike_timing_rcd{i,j} = 0.1*spike_timing_rcd{i,j};
        mean_rcd(i,j) = mean(ISI_rcd{i,j});
        cv_rcd(i,j) = std(ISI_rcd{i,j})/mean_rcd(i,j);
        y_rcd{i,j} = full(y_sparse_rcd{i,j});
        
        % GLM
        ISI = ISI_rcd{i,j};
        y = y_rcd{i,j};
        T = tot_t;
        I = zeros(T,1);
        y_glm = y(1,1:T);

        nkt = 100; % number of ms in stim filter
        kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
        kbasprs.ncos = 3; % number of raised-cosine vectors to use
        kbasprs.kpeaks = [5 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
        kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
        % basis functions for post-spike kernel
        ihbasprs.ncols = 10;  % number of basis vectors for post-spike kernel
        hPeaksMax = 3*max([mean(ISI) tau_E/10 tau_I/10 tau_M/10]);
        ihbasprs.hpeaks = [10 10*hPeaksMax];  % peak location for first and last vectors, in ms
        ihbasprs.b = 10*0.3*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
        ihbasprs.absref = 0; % absolute refractory period, in ms

        fit_k = 1;
        if max(I) == 0
            fit_k = 0;
        end
        plotFlag = 0;    % plot fit

        [k, h, dc, prs, kbasis, hbasis] = fit_glm(I,y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);
        history_rcd{i,j} = h;
        bias{i,j} = dc;
        
        % hazard
        ihbasprs.ncols = 15;  % number of basis vectors for post-spike kernel
        hPeaksMax = 10*max([mean(ISI) tau_E/10 tau_I/10 tau_M/10]);
        ihbasprs.hpeaks = [5 10*hPeaksMax];  % peak location for first and last vectors, in ms
        ihbasprs.b = 10*0.3*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
        ihbasprs.absref = 0; % absolute refractory period, in ms

        [k, h, dc, prs, kbasis, hbasis] = fit_hazard(I,y_glm',dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag);        
        hazard_rcd{i,j} = h;
        
    end
end
%% Plot filter & hazard function
figure
nn = 0;
for i = 1:length(tau_E_range)
    for j = 1:length(tau_M_range)
        nn = nn + 1;
        subplot(length(tau_E_range),length(tau_M_range),nn);
        plot_curve = history{i,j};
        plot((1:length(plot_curve))/10,plot_curve);
        xlim([1 length(plot_curve)/10 ]);
        title(['\tau_{input} = ',num2str(tau_E_range(i)),';\tau_M = ',num2str(tau_M_range(j))]);
        grid on
    end
end
figure
nn = 0;
for i = 1:length(tau_E_range)
    for j = 1:length(tau_M_range)
        nn = nn + 1;
        subplot(length(tau_E_range),length(tau_M_range),nn);
        plot_curve = 10*exp(hazard{i,j}(2:end));
        plot((1:length(plot_curve))/10,plot_curve);
        xlim([1 200]);
        title(['\tau_{input} = ',num2str(tau_E_range(i)),';\tau_M = ',num2str(tau_M_range(j))]);
        grid on
    end
end
figure
T = (10*3e3);
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
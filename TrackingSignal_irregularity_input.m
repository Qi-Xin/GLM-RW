clearvars;

tau_E = 1e-3;           % 1ms
tau_I = tau_E;
tau_M = 20;
dt = 1;
p = 5e-1;
q = 5e-1;

tot_t = 1e3;
bin = 1;   %ms
ddt = bin;
V_E = 0.08;
V_I = 0.08;
adjStepOrNot = 0;
adjValue = 50;
V_th = 1;
V_reset = 0;

maxSig = 1e-1;
signalType = 3; % 1 for no signal, 2 for square wave, 3 for inverse gamma
I_per = zeros(1,1e3);
if signalType == 2
    I_per(1:1000) = 1;
    I_per = I_per/max(I_per)*maxSig;
end
if signalType == 3
    pd = makedist('InverseGaussian','mu',2e2,'lambda',6e2);
    I_per(201:1000) = pdf(pd,[1:800]);
    I_per = I_per/max(I_per)*maxSig;
end

repnum = ceil(tot_t/length(I_per));
I = repmat(I_per,1,repnum);
I = I(1:tot_t);
I_noInp = zeros(1,tot_t);

%%
% Simulation
V_E = 0.08;
V_I = V_E;
N_input = 30;
mu = 20;
ratio = 0.8;
lambda = 1./ratio^2*mu;
cv = sqrt(lambda/mu);
y = zeros(1e3, 1e3);
isi_rec = [];
for i = 1:1e3
    [ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI_InvGammaInput(...
        N_input,tau_E,tau_I,tau_M,V_E,V_I,mu,lambda,V_th,V_reset,I,tot_t,dt);
    y(i,:) = full(y_sparse);
    isi_rec = [isi_rec, ISI];
end

% Plot Raster
T = 1e3;
N = 100;

figure
subplot(2,1,1);
plotraster(y(1:N,:),1:1e3,'Simulated Result');
title(['Spike train (num of input neuron=',num2str(N_input),') (input CV = ',num2str(ratio),')' ]);
xlabel('ms');
ylabel('trial');

% Tracking Signal
subplot(2,1,2);
hold on
tot_N = 1e3;
fr = sum(y);
if ddt ~= 1
    fr2 = reshape(fr,ddt,[]);
    fr2 = sum(fr2);
else
    fr2 = fr;
end
yyaxis left
semilogy( ((1:length(fr2))-1/2)*ddt ,fr2/tot_N/ddt );
% axis([0 5000 3 9]);
xlabel('t');
ylabel('Spike Probability');
yyaxis right
plot(repmat(I_per,1,5));
title('Spike Probability & Signal');
axis([0 1000 0 0.10]);
xlabel('t');
ylabel('Signal');

figure
edge = [0:10:5e2];
histogram(isi_rec,edge,'Normalization','probability');
ylim([0 0.5]);
title(['Histogram (num of input neuron=',num2str(N_input),') (input CV = ',num2str(ratio),')' ]);
xlabel('ISI (ms)');
ylabel('probability');

%% tracking accuracy of different CVs
ratio_list = [0.05,0.1,0.2,0.3,0.4:0.3:1.9] ;
ratio_list = [0.1,0.2,0.3,0.4:0.3:1.9] ;
n_rep = 5;
result = zeros(n_rep, length(ratio_list));
cv_output = zeros(1,length(ratio_list));
fr_base = zeros(1,length(ratio_list));
tot_N = 1e3;
% I = I+1e-2;
N_input = 1e3;
V_E = 0.02;
V_I = 0.02;


for i=1:length(ratio_list)
    
    ratio = ratio_list(i);
    mu = 20;
    lambda = 1./ratio^2*mu;
    for j=1:n_rep
        [i,j]
        cv = sqrt(mu/lambda);
        y = zeros(1e3, 1e3);
        y_sparse = {};
        parfor ii = 1:1e3
            [~,~,y_sparse{ii},~,~,~] = GetISI_InvGammaInput(...
                N_input,tau_E,tau_I,tau_M,V_E,V_I,mu,lambda,V_th,V_reset,I,tot_t,dt);
            y(ii,:) = full(y_sparse{ii});
        end
        fr = sum(y)/tot_N;
        p = polyfit(fr, I_per,1);
        yfit = polyval(p,fr);
        yresid = I_per - yfit;
        SS = sum(yresid.^2);
        result(j,i) = SS;
    end
    % sample ISIs during no stimulus
    isi_rec = [];
    epoch = 0;
    while true
        [ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI_InvGammaInput(...
            N_input,tau_E,tau_I,tau_M,V_E,V_I,mu,lambda,V_th,V_reset,I_noInp,tot_t,dt);
        isi_rec = [isi_rec, ISI];
        epoch = epoch + 1;
        if (length(isi_rec)>=1e4) | (epoch>=1000 & length(isi_rec)>=2e1)
            fr_base(i) = length(isi_rec)/epoch/1e3;
            break
        end
    end
    cv_output(i) = sqrt(var(isi_rec))/mean(isi_rec);
    length(isi_rec)
end

%% 
figure
subplot(3,1,1)
hold on
% yyaxis left
errorbar(ratio_list, mean(result), var(result)/sqrt(n_rep),'horizontal');

ylabel('Tracking ability');
title('Tracking accuracy (as sum of residual)');
legend(['num of input neuron=',num2str(N_input)]);

subplot(3,1,2)
% yyaxis right
plot(ratio_list, cv_output);
title('CV of output spike trains');
ylabel('CV');
legend(['num of input neuron=',num2str(N_input)]);

subplot(3,1,3)
% yyaxis right
plot(ratio_list, fr_base);
title('Baseline firing rate');
ylabel('FR');
legend(['num of input neuron=',num2str(N_input)]);
xlabel('CV of input spike trains');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% BELOW ARE POISSON %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
p = 300 *50/1e3;
q = p;
tot_N = 1e3;
V_E = 0.02;
V_I = V_E;

cv = 1;
y = zeros(1e3, 1e3);
y_sparse = {};
for ii = 1:tot_N
    [~,~,y_sparse{ii},~,~,~] = GetISI(...
        tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt);
    y(ii,:) = full(y_sparse{ii});
end
fr = sum(y)/tot_N;
poly = polyfit(fr, I_per,1);
yfit = polyval(poly,fr);
yresid = I_per - yfit;
SS = sum(yresid.^2);

% sample ISIs during no stimulus
isi_rec = [];
epoch = 0;
while true
    [ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI(...
        tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I_noInp,tot_t,dt);
    isi_rec = [isi_rec, ISI];
    epoch = epoch + 1;
    if (length(isi_rec)>=1e4) | (epoch>=1000 & length(isi_rec)>=2e1) | epoch>=1e4
        fr_base(i) = length(isi_rec)/epoch/1e3;
        break
    end
end
cv_output = sqrt(var(isi_rec))/mean(isi_rec)
length(isi_rec)

figure
subplot(2,1,1);
plotraster(y(1:N,:),1:1e3,'Simulated Result');
title(['Spike train (num of input neuron=',num2str(N_input),') (Poisson input, input CV = 1)' ]);
xlabel('ms');
ylabel('trial');

% Tracking Signal
subplot(2,1,2);
hold on
tot_N = 1e3;
fr = sum(y);
if ddt ~= 1
    fr2 = reshape(fr,ddt,[]);
    fr2 = sum(fr2);
else
    fr2 = fr;
end
yyaxis left
semilogy( ((1:length(fr2))-1/2)*ddt ,fr2/tot_N/ddt );
% axis([0 5000 3 9]);
xlabel('t');
ylabel('Spike Probability');
yyaxis right
plot(repmat(I_per,1,5));
title('Spike Probability & Signal');
axis([0 1000 0 0.10]);
xlabel('t');
ylabel('Signal');

%% tracking accuracy of different CVs under Poisson input
p_list = 50*[10:10:100,120:20:300]/1e3;
p_list = logspace(-1.3,4,20);
n_rep = 5;
result = zeros(n_rep, length(p_list));
cv_output = zeros(1,length(p_list));
fr_base = zeros(1,length(p_list));
tot_N = 1e3;
% I = I+1e-2;
N_input = 1e3;
V_E = 0.02;
V_I = V_E;

tot_t = 1e6;
repnum = ceil(tot_t/length(I_per));
I_multi = repmat(I_per,1,repnum);

for i=1:length(p_list)
    p = p_list(i);
    q = p;
    mu = 20;
    for j=1:n_rep
        [i,j]

        [ISI,spike_timing,y_sparse,V,inputE,inputI] = ... 
            GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I_multi,tot_t,dt);
        y = full(y_sparse);
        ST = (reshape(y,[],tot_N)');
        y = ST;
        
        fr = sum(y)/tot_N;
        poly = polyfit(fr, I_per,1);
        yfit = polyval(poly,fr);
        yresid = I_per - yfit;
        SS = sum(yresid.^2);
        result(j,i) = SS;
    end
    % sample ISIs during no stimulus
    isi_rec = [];
    epoch = 0;
    while true
        [ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI(...
            tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,repmat(I_noInp,1,repnum),tot_t,dt);
        isi_rec = [isi_rec, ISI];
        epoch = epoch + 1;
        if (length(isi_rec)>=1e4) | ((epoch>=1000 & length(isi_rec)>=2e1)) | epoch>=1e1
            fr_base(i) = length(isi_rec)/epoch/1e3;
            break
        end
    end
    cv_output(i) = sqrt(var(isi_rec))/mean(isi_rec);
    length(isi_rec)
end

%%
figure
subplot(3,1,1)
hold on
% yyaxis left
errorbar(p_list, mean(result), var(result)/sqrt(n_rep),'horizontal');

ylabel('Tracking ability');
title('Tracking accuracy (as sum of residual)');
legend(['num of input neuron=',num2str(N_input)]);

subplot(3,1,2)
% yyaxis right
plot(p_list, cv_output);
title('CV of output spike trains');
ylabel('CV');
legend(['num of input neuron=',num2str(N_input)]);

subplot(3,1,3)
% yyaxis right
plot(p_list, fr_base);
title('Baseline firing rate');
ylabel('FR');
legend(['num of input neuron=',num2str(N_input)]);
xlabel('CV of input spike trains');
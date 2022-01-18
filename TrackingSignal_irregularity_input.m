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

% Simulation
mu = 20;
lambda = 0.2*mu;
cv = sqrt(lambda/mu);
y = zeros(1e3, 1e3);
for i = 1:1e3
    [ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI_InvGammaInput(tau_E,tau_I,tau_M,V_E,V_I,mu,lambda,V_th,V_reset,I,tot_t,dt);
    y(i,:) = full(y_sparse);
end

% Plot Raster
T = 1e3;
N = 100;

figure
subplot(2,1,1);
plotraster(y(1:N,:),1:1e3,'Simulated Result');
title('Spike train');
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


%% tracking accuracy of different CVs
ratio_list = 0.1:0.3:1.9 ;
n_rep = 3;
result = zeros(n_rep, length(ratio_list));
cv_output = zeros(1,length(ratio_list));
for i=1:length(ratio_list)
    ratio = ratio_list(i);
    mu = 20;
    lambda = 1./ratio^2*mu;
    isi_rec = [];
    for j=1:n_rep
        [i,j]
        cv = sqrt(mu/lambda);
        y = zeros(1e3, 1e3);
        for ii = 1:1e3
            [ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI_InvGammaInput(tau_E,tau_I,tau_M,V_E,V_I,mu,lambda,V_th,V_reset,I,tot_t,dt);
            if ii<=10
                isi_rec = [isi_rec,ISI];
            end
            y(ii,:) = full(y_sparse);
        end
        fr = sum(y)/tot_N;
        p = polyfit(fr, I_per,1);
        yfit = polyval(p,fr);
        yresid = I_per - yfit;
        SS = sum(yresid.^2);
        result(j,i) = SS;
    end
    cv_output(i) = sqrt(var(isi_rec))/mean(isi_rec);
end
%% 
figure
hold on
yyaxis left
errorbar(ratio_list, mean(result), var(result)/sqrt(n_rep),'horizontal');
xlabel('CV of input spike trains');
ylabel('Tracking ability (as sum of residual)');
yyaxis right
plot(ratio_list, cv_output);
clearvars;

tau_E = 1e-3;           % 1ms
tau_I = tau_E;
tau_M = 20;
dt = 1;
p = 5e-1;
q = 5e-1;

tot_t = 5e3;
bin = 1;   %ms
ddt = bin;
V_E = 0.01/5;
V_I = 0.01/5;
adjStepOrNot = 0;
adjValue = 50;
V_th = 1;
V_reset = 0;

maxSig = 0.4e-1;
signalType = 3; % 1 for no signal, 2 for square wave, 3 for gamma
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
I = I(1:tot_t) + 5e-2;
I_noInp = zeros(1,tot_t);



% Simulation

y = zeros(1e4, 5e3);
for i = 1:1e4
    [ISI,spike_timing,y_sparse,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt);
    y(i,:) = full(y_sparse);
end


% Plot Raster
T = 1e3;
N = 100;

figure
subplot(2,1,1);
plotraster(y(1:N,:),1:5e3,'Simulated Result');
title('Spike train');
xlabel('ms');
ylabel('trial');


% Tracking Signal
subplot(2,1,2);
hold on
tot_N = 1e4;
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
axis([0 5000 0 0.10]);
xlabel('t');
ylabel('Signal');
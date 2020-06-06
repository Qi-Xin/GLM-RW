clearvars;

tau_E = 50;
tau_I = 50;
tau_M = 100;
dt = 1;
V_E = 1e-1*(1-exp(-dt/tau_E));
V_I = 1e-1*(1-exp(-dt/tau_I));
p = 2e-1;
q = 1e-1;

bin = 1;
V_th = 1;
V_reset = 0;
tot_t = 20*tau_M;
tot_neuron = 1e5;
V = NaN*zeros(tot_neuron,tot_t);
V(:,1) = V_reset; 
rdm = rand(tot_neuron,tot_t + 5*max(tau_E,tau_I));
inputE = zeros(tot_neuron,tot_t + 5*max(tau_E,tau_I));
inputI = zeros(tot_neuron,tot_t + 5*max(tau_E,tau_I));
tot_E = zeros(tot_neuron,1);
tot_I = zeros(tot_neuron,1);
%spike_timing = [0];

for t = 2: ( tot_t + 5*max(tau_E,tau_I) )
    inputE(:,t) = inputE(:,t-1).*exp(-dt/tau_E);
    inputI(:,t) = inputI(:,t-1).*exp(-dt/tau_I);
    inputE(:,t) = inputE(:,t) + (rdm(:,t) <= p) .* V_E;
    inputI(:,t) = inputI(:,t) + (rdm(:,t) >= 1-q) .* V_I;
end
inputE = inputE(: , (5*max(tau_E,tau_I)+1) :end);
inputI = inputI(: , (5*max(tau_E,tau_I)+1) :end);
for t = 2:tot_t
    V(:,t) = V(:,t-1).*exp(-dt/tau_M);
    V(:,t) = V(:,t-1) + inputE(:,t) - inputI(:,t);
    temp = V(:,t);
    temp(find(temp<0)) = 0;
    V(:,t) = temp;
end

%ISI = diff(spike_timing);

%%
plotT = [0.2*tau_E,0.4*tau_E,0.6*tau_E,1*tau_E,2*tau_E,3*tau_E];
Vgrid = 0:(5*V_E):(1e3*V_E);
rcd = NaN*ones(length(plotT),length(Vgrid)-1);

figure
hold on
for i = 1:length(plotT)
    h = histogram(V(:,plotT(i)),Vgrid,'Normalization','pdf','DisplayStyle','stairs');
    rcd(i,:) = h.Values;
end

Vgrid = Vgrid(1:end-1);
figure
hold on
for i = 1:length(plotT)
    plot(Vgrid,rcd(i,:));
end
title('P(x,t=t0)');
legend('t=0.1*tau','t=0.5*tau','t=1*tau','t=2*tau','t=5*tau','t=10*tau');
xlabel('steps');
ylabel('Probability');
axis([0 5 0 5]);


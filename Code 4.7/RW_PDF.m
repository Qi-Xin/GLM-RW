clearvars;

tau_M = 100;
V_E = 1;
V_I = 1;
p = 2e-1;
q = 1e-1;

bin = 1;
V_th = 0;
V_reset = 1;
tot_t = 20*tau_M;
tot_neuron = 1e5;
V = NaN*zeros(tot_neuron,tot_t);
V(:,1) = V_reset; 
rdm = rand(tot_neuron,tot_t);
tot_E = zeros(tot_neuron,1);
tot_I = zeros(tot_neuron,1);
dt = 1;
%spike_timing = [0];

for t = 2:tot_t
    V(:,t) = V(:,t-1).*exp(-dt/tau_M);
    V(:,t) = V(:,t) + (rdm(:,t) <= p).*V_E - (rdm(:,t) >= 1-q).*V_I;
    temp = V(:,t);
    temp(find(temp<0)) = 0;
    V(:,t) = temp;
end

%ISI = diff(spike_timing);

%%
plotT = [1,0.1*tau_M,0.5*tau_M,tau_M,2*tau_M,5*tau_M,10*tau_M];
Vgrid = 0:V_E:(7e2*V_E);
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
    plot(Vgrid-1,rcd(i,:));
end
title('P(x,t=t0)');
legend('t=0','t=0.1*tau','t=0.5*tau','t=1*tau','t=2*tau','t=5*tau','t=10*tau');
xlabel('steps');
ylabel('Probability');
axis([-2 50 0 1]);


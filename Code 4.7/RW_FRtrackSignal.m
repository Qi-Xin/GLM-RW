clearvars;
rng(3)

tau_M = 100;
tau_E = 100;
tau_I = 100;
p = 2e-1;
q = 1e-1;
V_E = 2e-1;
V_I = 1e-1;
% for balanced / low / high noise
%{
p = 2e-1;
q = 1e-1;
V_E = 2e-1;
V_I = 2e-1;

p = 2e-1/7;
q = 1e-1/7;
V_E = 2e-1;
V_I = 2e-1;

p = 2e-1*3;
q = 1e-1*3;
V_E = 2e-1;
V_I = 2e-1;
%}
tot_t = 5e3;
I = zeros(1,tot_t);
for t = 2:tot_t
    I(t) = I(t-1)+1e-2*(rand(1,1)-0.5);
end
I = smooth(I,100);
I = 0.03*I;
ddt = 10;

bin = 10;
V_th = 1;
V_reset = 0;
tot_neuron = 1e4;
V = NaN*zeros(tot_neuron,tot_t);
V(:,1) = V_reset; 
rdm = rand(tot_neuron,tot_t);
tot_E = zeros(tot_neuron,1);
tot_I = zeros(tot_neuron,1);
dt = 1;
spike_timing = zeros(tot_neuron,tot_t);

for t = 2:tot_t
    V(:,t) = V(:,t-1).*exp(-dt/tau_M);
    V(:,t) = V(:,t) + (rdm(:,t) <= p).*V_E - (rdm(:,t) >= 1-q).*V_I + I(t);
    spike_timing(:,t) = ( V(:,t)>=1 );
    temp = V(:,t);
    temp(find(temp<0)) = 0;
    V(:,t) = temp;
    temp = V(:,t);
    temp(find(temp>1)) = 0;
    V(:,t) = temp;
end

%ISI = diff(spike_timing);

%%
figure
hold on
fr = sum(spike_timing);
fr2 = reshape(fr,[ddt,tot_t/ddt]);
fr2 = sum(fr2);
yyaxis left
semilogy( ((1:length(fr2))-1/2)*ddt ,fr2/tot_neuron*1e2);
%axis([0 5000 3 9]);
%axis([0 5000 55 65]);
xlabel('t');
ylabel('Firing rate');
yyaxis right
plot(1:tot_t,I);
title('Firing rate & Signal');
%axis([0 5000 -0.15 0.25]);
xlabel('t');
ylabel('Signal');


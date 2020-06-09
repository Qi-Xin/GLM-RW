function [ISI,spike_timing,y,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt)

V = NaN*zeros(1,tot_t);
V(1) = 0; 
y = zeros(1,tot_t);
spike_timing = [0];

LB = -0.1;
ratio = 0.8;
tau_E1 = 5;
tau_E2 = 67;
tau_I1 = 5;
tau_I2 = 150;

add = ceil(5*max([tau_E,tau_I,1e3/dt]));
rdmE1 = random('poisson',p*ratio,1,tot_t+add);
rdmE2 = random('poisson',p*(1-ratio),1,tot_t+add);
rdmI1 = random('poisson',q*ratio,1,tot_t+add);
rdmI2 = random('poisson',q*(1-ratio),1,tot_t+add);
inputE1 = zeros(1,tot_t + add);
inputE2 = zeros(1,tot_t + add);
inputI1 = zeros(1,tot_t + add);
inputI2 = zeros(1,tot_t + add);

for t = 2:( tot_t + add )
	inputE1(t) = inputE1(t-1)*exp(-dt/tau_E1);
    inputI1(t) = inputI1(t-1)*exp(-dt/tau_I1);
    inputE2(t) = inputE2(t-1)*exp(-dt/tau_E2);
    inputI2(t) = inputI2(t-1)*exp(-dt/tau_I2);
    inputE1(t) = inputE1(t) + V_E*rdmE1(t);
    inputI1(t) = inputI1(t) + V_I*rdmI1(t);
    inputE2(t) = inputE2(t) + V_E*rdmE2(t);
    inputI2(t) = inputI2(t) + V_I*rdmI2(t);
end
inputE1 = inputE1( (add+1) :end);
inputI1 = inputI1( (add+1) :end);
inputE2 = inputE2( (add+1) :end);
inputI2 = inputI2( (add+1) :end);

for t = 2:tot_t
    if V(t-1)>=1
        spike_timing = [spike_timing,t-1];
        y(t-1) = 1;
        V(t) = 0;
        continue
    end
    
    V(t) = V(t-1)*exp(-dt/tau_M);
    V(t) = V(t) + inputE1(t) - inputI1(t) + inputE2(t) - inputI2(t) + I(t);
    
    if V(t)<= LB
        V(t) = LB;
    end
    
end

ISI = diff(spike_timing);
y = sparse(y);
inputE = inputE1+inputE2;
inputI = inputI1+inputI2;
end


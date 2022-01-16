function [ISI,spike_timing,y,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt,LB)

a = p;
b = q;

V = NaN*zeros(1,tot_t);
V(1) = 0; 
y = zeros(1,tot_t);
spike_timing = [0];

add = ceil(5*max([tau_E,tau_I,1e3/dt]));

rdmE = zeros(1,tot_t + add);
t = 0;
while t<=tot_t+add
    deltat = 1/gamrnd(a,b);
    t = t + deltat;
    rdmE(ceil(t)) = 1;
end

rdmI = zeros(1,tot_t + add);
t = 0;
while t<=tot_t+add
    deltat = 1/gamrnd(a,b);
    t = t + deltat;
    rdmI(ceil(t)) = 1;
end

inputE = zeros(1,tot_t + add);
inputI = zeros(1,tot_t + add);

if ~exist('LB')
    LB = -V_E;
end
for t = 2:( tot_t + add )
	inputE(t) = inputE(t-1)*exp(-dt/tau_E);
    inputI(t) = inputI(t-1)*exp(-dt/tau_I);
    inputE(t) = inputE(t) + V_E*rdmE(t);
    inputI(t) = inputI(t) + V_I*rdmI(t);
end
inputE = inputE( (add+1) :end);
inputI = inputI( (add+1) :end);
sum(inputE)
for t = 2:tot_t
    if V(t-1)>=V_th
        spike_timing = [spike_timing,t-1];
        y(t-1) = 1;
        V(t) = V_reset;
        continue
    end
    
    V(t) = V(t-1)*exp(-dt/tau_M);
    V(t) = V(t) + inputE(t) - inputI(t) + I(t);
    
    if V(t)<= LB
        V(t) = LB;
    end
    
end

ISI = diff(spike_timing);
y = sparse(y);

end

function [ISI,spike_timing,y,V,inputE,inputI] = GetISI_VaryingNoise(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt,LB)

V = NaN*zeros(1,tot_t);
V(1) = 0; 
y = zeros(1,tot_t);
spike_timing = [];

add = ceil(5*max([tau_E,tau_I,1e3/dt]));
rdmE = random('poisson',p,1,tot_t+add);
rdmI = random('poisson',q,1,tot_t+add);
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

v = 0;
for t = 1:add
    if v>=V_th
        v = V_reset;
        continue
    end
    
    v = v*exp(-dt/tau_M);
    v = v + inputE(t) - inputI(t);
    
    if v<= LB
        v = LB;
    end
    
end
inputE = inputE( (add+1) :end);
inputI = inputI( (add+1) :end);
V(1) = v;

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

if V(t)>=V_th
    spike_timing = [spike_timing,t];
    y(t) = 1;
end

ISI = diff(spike_timing);
y = sparse(y);

end


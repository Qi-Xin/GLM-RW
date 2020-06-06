function [ISI,spike_timing,y,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,tot_N,tot_t,dt)

V = NaN*zeros(1,tot_t);
V(1) = 0; 
y = zeros(1,tot_t);
spike_timing = [0];

add = ceil(5*max([tau_E,tau_I,1e3/dt]));
rdm = rand(1,tot_t+add);
inputE = zeros(1,tot_t + add);
inputI = zeros(1,tot_t + add);
for t = 2:( tot_t + add )
	inputE(t) = inputE(t-1)*exp(-dt/tau_E);
    inputI(t) = inputI(t-1)*exp(-dt/tau_I);
    if rdm(t) <= p
        inputE(t) = inputE(t) + V_E;
    elseif rdm(t) >= 1-q
        inputI(t) = inputI(t) + V_I;
    end
end
inputE = inputE( (add+1) :end);
inputI = inputI( (add+1) :end);
for t = 2:tot_t
    V(t) = V(t-1)*exp(-dt/tau_M);
    V(t) = V(t) + inputE(t) - inputI(t);

    if V(t)>=1
        spike_timing = [spike_timing,t];
        y(t) = 1;
        if length(spike_timing)>=tot_N
            break;
        end
        V(t) = 0;
    end

    if V(t)<=0
        V(t) = 0;
    end

end

ISI = diff(spike_timing);
y = sparse(y);

end


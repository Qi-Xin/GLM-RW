function [ISI,spike_timing,y,V,inputE,inputI] = GetISI(tau_E1,tau_E2,tau_I1,tau_I2,ratio,tau_M,V_E,V_I,p,q,V_th,V_reset,tot_N,tot_t,dt)

V = NaN*zeros(1,tot_t);
V(1) = 0; 
y = zeros(1,tot_t);
spike_timing = [0];

p1 = p*(ratio);
p2 = p*(1-ratio);
q1 = q*(ratio);
q2 = q*(1-ratio);
%add = ceil(5*max(tau_E,tau_I));
add = 0;
rdm1 = rand(1,tot_t+add);
rdm2 = rand(1,tot_t+add);
inputE1 = zeros(1,tot_t + add);
inputE2 = zeros(1,tot_t + add);
inputI1 = zeros(1,tot_t + add);
inputI2 = zeros(1,tot_t + add);
for t = 2:( tot_t + add )
	inputE1(t) = inputE1(t-1)*exp(-dt/tau_E1);
    inputE2(t) = inputE2(t-1)*exp(-dt/tau_E2);
    inputI1(t) = inputI1(t-1)*exp(-dt/tau_I1);
    inputI2(t) = inputI2(t-1)*exp(-dt/tau_I2);
    if rdm1(t) <= p1
        inputE1(t) = inputE1(t) + V_E;
    elseif rdm1(t) >= 1-q1
        inputI1(t) = inputI1(t) + V_I;
    end
    if rdm2(t) <= p2
        inputE2(t) = inputE2(t) + V_E;
    elseif rdm2(t) >= 1-q2
        inputI2(t) = inputI2(t) + V_I;
    end
end
inputE1 = inputE1( (add+1) :end);
inputE2 = inputE2( (add+1) :end);
inputI1 = inputI1( (add+1) :end);
inputI2 = inputI2( (add+1) :end);
for t = 2:tot_t
    V(t) = V(t-1)*exp(-dt/tau_M);
    V(t) = V(t) + inputE1(t) - inputI1(t) + inputE2(t) - inputI2(t);

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


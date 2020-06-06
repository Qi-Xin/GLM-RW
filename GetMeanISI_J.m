function [J] = GetMeanISI_J(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt,adjValue)

[ISI,spike_timing] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt);
upper = (mean(ISI)+std(ISI)/sqrt(length(ISI)));
lower = (mean(ISI)-std(ISI)/sqrt(length(ISI)));
J = mean(ISI) - adjValue;

if length(spike_timing) == 1
    J = Inf;
end

end

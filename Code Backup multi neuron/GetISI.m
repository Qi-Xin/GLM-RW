function [ISI,spike_timing,y,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,tot_N,tot_t,dt,Nn)

V = NaN*zeros(Nn,tot_t);
V(:,1) = zeros(Nn,1); 
y = zeros(Nn,tot_t);
spike_timing = zeros(Nn,1);
spike_number = zeros(Nn,1);

add = ceil(10*max([tau_E,tau_I]));
rdm = rand(Nn,tot_t+add);
inputE = zeros(Nn,tot_t + add);
inputI = zeros(Nn,tot_t + add);

for t = 2:( tot_t + add )
	inputE(:,t) = inputE(:,t-1)*exp(-dt/tau_E);
    inputI(:,t) = inputI(:,t-1)*exp(-dt/tau_I);
    inputE(:,t) = inputE(:,t) + V_E*(rdm(:,t)<=p);
    inputI(:,t) = inputI(:,t) + V_I*(rdm(:,t)>=1-q);
end

inputE = inputE( : , (add+1):end);
inputI = inputI( : , (add+1):end);
for t = 2:tot_t
    V(:,t) = V(:,t-1)*exp(-dt/tau_M);
    V(:,t) = V(:,t) + inputE(:,t) - inputI(:,t);
    
    havesp = find(V(:,t)>=1);
    if length(havesp) ~= 0
        
        for j = 1:length(havesp)
            Nnum = havesp(j);
            spike_number(Nnum) = spike_number(Nnum)+1;
            spike_timing(Nnum,spike_number(Nnum)) = t;
            y(Nnum,t) = 1;
            V(Nnum,t) = 0;
        end
    end
    
    V(:,t) = ( abs(V(:,t)) + V(:,t) )/2;

end

ISI = diff(spike_timing,1,2);
ISI = ISI(find(ISI>0));
y = sparse(y);

end


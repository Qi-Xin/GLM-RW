clearvars;

tau_E = 500;
tau_I = tau_E;
tau_M = 100;
dt = 1;
V_E = 1e-1*(1-exp(-dt/tau_E));
V_I = 1e-1*(1-exp(-dt/tau_I));
p = 2e-1;
q = 1e-1;

bin = 20;
V_th = 1;
V_reset = 0;
tot_t = 1e6;
V = NaN*zeros(1,tot_t);
V(1) = 0; 
spike_timing = [0];

%add = ceil(5*max(tau_E,tau_I));
add = 0;
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
        V(t) = 0;
    end
   
    if V(t)<=0
        V(t) = 0;
    end
    
end

ISI = diff(spike_timing);

%%
ddt = bin;
max1 = ceil(max(ISI));
dist = zeros(1,ceil(max1/ddt));
total_trial = length(ISI);
for j = 1:total_trial
    dist(ceil(ISI(j)/ddt)) = dist(ceil(ISI(j)/ddt))+1;
end

dead = 0;
svv(1) = total_trial;
for i = 1:length(dist)
    dead = dead + dist(i);
    svv(i+1) = total_trial-dead;
end
svv = svv(1:length(svv)-1);
[rate,err] = binofit(dist,svv);
xx1 = ((1:length(dist))-1/2);
xdata = (xx1)*ddt;
ydata = rate;
last = max(find(diff(err')<=0.1));
xdata = xdata(1:last);
ydata = ydata(1:last);

figure
yyaxis left
h = histogram(ISI,0:ddt:max1,'Normalization','pdf');
xlabel('t/ms');
ylabel('ISI distribution');
axis([0 1000 0 0.0075]);
yyaxis right
errorbar(xdata,rate(1:last)/ddt,err(1:last,1)'/ddt-rate(1:last)/ddt,-err(1:last,2)'/ddt+rate(1:last)/ddt);
ylabel('Hazard function');
axis([0 1000 0 0.015]);
title(['tau_E = ',num2str(tau_E), ...
    ';tau_I = ',num2str(tau_I), ...
    ';tau_M = ',num2str(tau_M), ...
    ';V_E = ',num2str(V_E), ...
    ';V_I = ',num2str(V_I), ...
    ';p = ',num2str(p), ...
    ';q = ',num2str(q)]);

%{
title('Hazard function');
legend('r2');
xlabel('t/s');
ylabel('rate');
axis([0 200 0 1]);
%}


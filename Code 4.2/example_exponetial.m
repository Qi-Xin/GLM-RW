clearvars;

tau_E = 20;
tau_I = 20;
V_E = 2e-1;
V_I = 2e-1;
p = 2e-1;
q = 1e-1;

bin = 20;
V_th = 0;
V_reset = 1;
tot_t = 1e7;
V = NaN*zeros(1,tot_t);
V(1) = V_reset; 
rdm = rand(1,tot_t);
tot_E = 0;
tot_I = 0;
dt = 1;
spike_timing = [0];

for t = 2:tot_t
    tot_E = tot_E*exp(-dt/tau_E);
    tot_I = tot_I*exp(-dt/tau_I);
    if rdm(t) <= p
        tot_E = tot_E + V_E;
    elseif rdm(t) >= 1-q
        tot_I = tot_I + V_I;
    end
    V(t) = tot_E - tot_I;
    if V(t)>=1
        spike_timing = [spike_timing,t];
        tot_E = 0;
        tot_I = 0;
    end
    if V(t)<=0
        tot_E = 0;
        tot_I = 0;
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
yyaxis right
errorbar(xdata,rate(1:last)/ddt,err(1:last,1)'/ddt-rate(1:last)/ddt,-err(1:last,2)'/ddt+rate(1:last)/ddt);
ylabel('Hazard function');
%axis([0 max1 0 1]);
title(['tau_E = ',num2str(tau_E), ...
    ';tau_I = ',num2str(tau_I), ...
    ';V_E = ',num2str(V_E), ...
    ';V_I = ',num2str(V_I), ...
    ';p = ',num2str(p), ...
    ';q = ',num2str(q)]);
axis([0 1000 0 0.02]);
%{
title('Hazard function');
legend('r2');
xlabel('t/s');
ylabel('rate');
axis([0 200 0 1]);
%}


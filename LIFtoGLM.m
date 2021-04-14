clearvars;

p_list = [1e-1,3e-1,1e0,2.5e0,0.5e1,1e1];
I_list = linspace(0,4.8e-2,20);
BalaV = NaN*zeros(length(p_list),length(I_list));
lambda = NaN*zeros(length(p_list),length(I_list));

tot_t = 1e6;
bin = 1;   %ms
ddt = bin;
V_E = 0.1;
V_I = 0.1;
V_th = 1;
V_reset = 0;
tau_E = 1e1;       % 1ms
tau_I = 1e1;
tau_M = 20;
dt = 1;

for ii = 1:length(p_list)
    for jj = 1:length(I_list)
        [ii,jj]
        p = p_list(ii);
        q = p;
        I_per = ones(1,1e3)*I_list(jj);
        repnum = ceil(tot_t/length(I_per));
        I = repmat(I_per,1,repnum);
        I = I(1:tot_t);
        
        % Get balanced position
        [~,~,~,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,0,0,V_th,V_reset,I_per,1e3,dt,-0.1);
        BalaV(ii,jj) = V(end);
        
        % firing rate
        [ISI,~,~,V,~,~] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,BalaV(ii,jj),I,tot_t,dt,-0.1);
        lambda(ii,jj) = length(ISI)/(tot_t/1e3);
    end
end
%%
figure;
hold on
for ii = 1:length(p_list)
    plot(BalaV(ii,:),log(lambda(ii,:)./1e3));
end
xlabel('Mean Voltage');
ylabel('Firing rate');
title({'Mean Integrated Squared Error';'between True Stimulus and Normalized PSTH'});

%%
p = 1e-3;
q = 0;
tot_t = 1e7;
bin = 1;
ddt = bin;
V_E = 0.1;
V_I = 0;
V_th = 1e8;
V_reset = 0;
tau_E = 1e-5;
tau_I = 1e-5;
tau_M = 1e2;
dt = 1;
[~,~,~,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,1e8,V_reset,0*ones(1,tot_t),tot_t,dt,-1);
figure
h = histogram(V,5e1);
yy = h.Values(10:end);
xx = h.BinEdges(11:end);

%%
T = 1e5;
V = ones(1,T);
V(1) = 1;
tau = 1e3;
for t = 2:T
    V(t) = V(t-1)*exp(-1/tau);
    if mod(t,1e4) == 0
        V(t) = V(t)+1;
    end
end
figure
h = histogram(V,5e1);
yy = h.Values(10:end-3);
xx = h.BinEdges(11:end-3);
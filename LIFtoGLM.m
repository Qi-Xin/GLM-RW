clearvars;

p_list = [1e0,2.5e0,0.5e1,1e1,1.5e1,2.5e1,5e1,15e1,3e2];
I_list = linspace(0,4.8e-2,20);
BalaV = NaN*zeros(length(p_list),length(I_list));
lambda = NaN*zeros(length(p_list),length(I_list));

tot_t = 1e6;
bin = 1;   %ms
ddt = bin;
V_E = 0.02;
V_I = 0.02;
V_th = 1;
V_reset = 0;
tau_E = 1e-3;       % 1ms
tau_I = tau_E;
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
        [~,~,~,V,~,~] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,0,0,V_th,V_reset,I_per,1e3,dt,-0.1);
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

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

%% Exactly happen to be exponential 
p = 5e-3;
q = 0;
tot_t = 1e8;
bin = 1;
ddt = bin;
V_E = 2e-1;
V_I = 0;
V_th = 1e8;
V_reset = 0;
tau_E = 1e-3;
tau_I = 1e-8;
tau_M = 1e2;
dt = 1;
[~,~,~,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,1e8,V_reset,0*ones(1,tot_t),tot_t,dt,-1);
figure
h = histogram(V,5e1,'Normalization','pdf');
yy = h.Values;
xx = h.BinEdges;
%xx = cumsum(xx);
%yy = cumsum(yy);
start = 25;
endnum = 0;
xx= xx(start+1:end-endnum)';
yy = yy(start:end-endnum)';
[f1,gof1] = fit(xx(1:end),yy(1:end),'exp1');
%[f2,gof2] = fit(xx(1:3),yy(1:3),'gauss1');
ft = fittype('b/(c*x+a)','independent','x','dependent','y');
[f3,gof3] = fit(xx,yy,ft);

figure
% subplot(2,2,1)
plot(f1,xx,yy);
title('Exp tail');
% subplot(2,2,2)
% plot(f2,xx,yy);
% title('Gaussian tail');
% subplot(2,2,3)
% plot(f3,xx,yy);
% title('1/x tail');

%% 1/x
p = 0;
q = 0;
tot_t = 1e8;
bin = 1;
ddt = bin;
V_E = 2e0;
V_I = 0;
V_th = 1e8;
I = zeros(1,tot_t);
I(1:5e4:end) = V_E;
V_reset = 0;
tau_E = 1e-3;
tau_I = 0e-8;
tau_M = 1e2;
dt = 0.1;
[~,~,~,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,1e8,V_reset,I,tot_t,dt,-1);
figure
h = histogram(V,5e1,'Normalization','pdf');
yy = h.Values;
xx = h.BinEdges;
%xx = cumsum(xx);
%yy = cumsum(yy);
start = 5;
endnum = 0;
xx= xx(start+1:end-endnum)';
yy = yy(start:end-endnum)';
[f1,gof1] = fit(xx(1:end),yy(1:end),'exp1');
%[f2,gof2] = fit(xx(1:3),yy(1:3),'gauss1');
ft = fittype('b/(c*x+a)','independent','x','dependent','y');
[f3,gof3] = fit(xx,yy,ft);

figure
% subplot(2,2,1)
plot(f1,xx,yy);
title('Exp tail');
% subplot(2,2,2)
% plot(f2,xx,yy);
% title('Gaussian tail');
% subplot(2,2,3)
plot(f3,xx,yy);
title('1/x tail');

figure;
plot(V(1:2e5));


%% Gaussian
p = 5e-2;
q = 0;
tot_t = 1e8;
bin = 1;
ddt = bin;
V_E = 1e-1;
V_I = 0;
V_th = 1e8;
V_reset = 0;
tau_E = 1e-3;
tau_I = 1e-8;
tau_M = 1e2;
dt = 1;
[~,~,~,V,inputE,inputI] = GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,1e8,V_reset,0*ones(1,tot_t),tot_t,dt,-1);
figure
h = histogram(V,5e1,'Normalization','pdf');
yy = h.Values;
xx = h.BinEdges;
%xx = cumsum(xx);
%yy = cumsum(yy);
start = 30;
endnum = 0;
xx= xx(start+1:end-endnum)';
yy = yy(start:end-endnum)';
[f1,gof1] = fit(xx(1:end),yy(1:end),'exp1');
[f2,gof2] = fit(xx(1:3),yy(1:3),'gauss1');
ft = fittype('b/(c*x+a)','independent','x','dependent','y');
[f3,gof3] = fit(xx,yy,ft);

figure
hold on
% subplot(2,2,1)
plot(f1,xx,yy);
title('Exp tail');
% subplot(2,2,2)
plot(f2,xx,yy);
title('Gaussian tail');
% subplot(2,2,3)
% plot(f3,xx,yy);
% title('1/x tail');
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
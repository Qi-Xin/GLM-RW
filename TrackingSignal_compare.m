%%% change line 7 p_list for range
%%% change line 44 for E/I balance or only E,
%%% change line 14 and 15 for steps

clearvars;
p_list = logspace(-1.3,4,20);     % balance
%p_list = logspace(-1.3,2,35);    % only Excitation
repeat = 10;
IMSE_raw = NaN*zeros(length(p_list),repeat);

tot_t = 1e6;
bin = 1;   %ms
ddt = bin;
V_E = 0.02;
V_I = V_E;
V_th = 1;
V_reset = 0;
tau_E = 1e-3;       % 1ms
tau_I = tau_E;
tau_M = 20;
dt = 1;

maxSig = 0.9e-1;
signalType = 3; % 1 for no signal, 2 for square wave, 3 for gamma
I_per = zeros(1,1e3);
if signalType == 2
    I_per(1:1000) = 1;
    I_per = I_per/max(I_per)*maxSig;
end
if signalType == 3
    pd = makedist('InverseGaussian','mu',2e2,'lambda',6e2);
    I_per(201:1000) = pdf(pd,[1:800]);
    I_per = I_per/max(I_per)*maxSig;
end
%plot(I_per);
repnum = ceil(tot_t/length(I_per));
I = repmat(I_per,1,repnum);
I = I(1:tot_t);
I_noInp = zeros(1,tot_t);

for ii = 1:length(p_list)
    p = p_list(ii);
    q = 0;
    %q = p;
    for jj = 1:repeat
        [ii,jj]
        % Simulation
        [ISI,spike_timing,y_sparse,V,inputE,inputI] = ... 
            GetISI(tau_E,tau_I,tau_M,V_E,V_I,p,q,V_th,V_reset,I,tot_t,dt);
        y = full(y_sparse);

        T = 1e3;
        N = 100;
        tot_N = (tot_t/length(I_per));
        fr = sum(reshape(y,[],tot_N)');
        if ddt ~= 1
            fr2 = reshape(fr,ddt,[]);
            fr2 = sum(fr2);
        else
            fr2 = fr;
        end

        % MISE
        target = fr2/tot_N/ddt;
        groundTrue = I_per;
        %baseline = mean(target(1:200));
        %fun = @(k) sum(((target-baseline)*k - groundTrue).^2)
        %x0 = 1;
        fun = @(k) sum(((target-k(2))*k(1) - groundTrue).^2)
        x0 = [1,0.2];
        [x,fval] = fminunc(fun,x0);
        IMSE_raw(ii,jj) = fval;
    end
end
%%
figure;
hold on
errorbar(p_list*1e3/50,mean(IMSE_raw,2),std(IMSE_raw')');
xlabel('Number of Background Noise Neurons (assuming each neuron is 50 Hz)');
ylabel('MISE');
title({'Mean Integrated Squared Error';'between True Stimulus and Normalized PSTH'});


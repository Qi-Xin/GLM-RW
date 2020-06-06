clearvars
%% 
AM = 2e-2;
tau_E = 200;
tau_I = 200;
V_E = 1e-1;
V_I = 1e-1;
p = 2e-1;
q = 1e-1;
tot_t = 5e3;
I = zeros(1,tot_t);
ddt = 10;

bin = 100;
V_th = 0;
V_reset = 1;
tot_neuron = 50;
V = NaN*zeros(tot_neuron,tot_t);
V(:,1) = V_reset; 
rdm = rand(tot_neuron,tot_t);
tot_E = zeros(tot_neuron,1);
tot_I = zeros(tot_neuron,1);
other = zeros(tot_neuron,1);
dt = 1;
spike_timing = zeros(tot_neuron,tot_t);

for t = 2:tot_t
    tot_E = tot_E.*exp(-dt/tau_E);
    tot_I = tot_I.*exp(-dt/tau_I);
    other = other.*exp(-dt/tau_E);
    tot_E = tot_E + (rdm(:,t) <= p).*V_E;
    other = other + (1-mod(ceil(t/1000),2))*AM;
    tot_I = tot_I + + (rdm(:,t) >= 1-q).*V_I;
    V(:,t) = tot_E - tot_I + other;
    spike_timing(:,t) = ( V(:,t)>=1 );
    tot_E(find( V(:,t)>=1 )) = 0;
    tot_I(find( V(:,t)>=1 )) = 0;
    other(find( V(:,t)>=1 )) = 0;
end

%ISI = diff(spike_timing);
for i = 1:tot_neuron
    pos = find(spike_timing(i,:) == 1);
    totalspikes(i) = length(pos);
    for j = 1:length(pos)
        sp(i,j)=pos(j);
    end
end

plotraster(spike_timing,1:tot_t,'Raster');

%%

figure
hold on
fr = sum(spike_timing);
fr2 = reshape(fr,[ddt,tot_t/ddt]);
fr2 = sum(fr2);
subplot(2,1,1)
plot( ((1:length(fr2))-1/2)*ddt ,fr2);
xlabel('t');
ylabel('Firing rate');
subplot(2,1,2)
plot(1:tot_t,(1-mod(ceil((1:tot_t)/1000),2))*AM);
axis([0 5000 0 5e-3]);
xlabel('t');
ylabel('Signal');

%%
plotornot=[0,1,1,0];     % whether to plot bases / the two filters / KS result / simulations 
kswhich  = 1;     %  1 for my filter, 0 for the original filter

nBases=5;
Baseslength=100;
tstrech=5;

Ntrial = tot_neuron;
T = tot_t;
SimulationT = 1e3;

trainingtrials=[1:tot_neuron];%1:10;

unitOfTime = 'ms';
binSize = 1; % 1 ms

%% Fit a model
expt = buildGLM.initExperiment(unitOfTime, binSize);
expt = buildGLM.registerSpikeTrain(expt, 'sptrain', 'simulated neuron');
expt2 = buildGLM.initExperiment(unitOfTime, binSize);
expt2 = buildGLM.registerSpikeTrain(expt2, 'sptrain', 'simulated neuron');
%% get y
for i=1:Ntrial
    expt.trial(i).sptrain = sp(i,1:totalspikes(i));
    expt.trial(i).duration = T;
end

dspec = buildGLM.initDesignSpec(expt);
bs = basisFactory.makeNonlinearRaisedCos(nBases, 1, [0 Baseslength], tstrech);

dspec = buildGLM.addCovariateSpiketrain(dspec, 'hist', 'sptrain', 'History filter', bs);

dm = buildGLM.compileSparseDesignMatrix(dspec, trainingtrials);
dm = buildGLM.addBiasColumn(dm);

filter_length=length(bs.tr);

y = buildGLM.getBinnedSpikeTrain(expt, 'sptrain', dm.trialIndices);


%% get dm2.X of input
for i=1:Ntrial
    expt2.trial(i).sptrain = [1001:2000,3001:4000]';
    expt2.trial(i).duration = T;
end

dspec2 = buildGLM.initDesignSpec(expt2);
bs2 = basisFactory.makeNonlinearRaisedCos(nBases, 1, [0 Baseslength], tstrech);

dspec2 = buildGLM.addCovariateSpiketrain(dspec2, 'hist', 'sptrain', 'History filter', bs2);

dm2 = buildGLM.compileSparseDesignMatrix(dspec2, trainingtrials);
dm2 = buildGLM.addBiasColumn(dm2);

filter_length=length(bs2.tr);


%% Do the regression
addpath('matRegress')
temp = sparse([full(dm2.X),full(dm.X)]);
wInit = temp \ y;
fnlin = @nlfuns.exp; % inverse link function (a.k.a. nonlinearity)
lfunc = @(w)(glms.neglog.poisson(w, temp, y, fnlin)); % cost/loss function

opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');

[wML, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
wML_full_both=full(wML);
%wvar = diag(inv(hessian));

%ws = buildGLM.combineWeights(dm2, wML);

%% Plot Bases
if plotornot(1)==1
    figure
    for i = 1:nBases
        plot(bs.B(:,i));
        hold on
    end
xlabel('t/ms','FontSize',13);
%ylabel();
title('history filter basis','FontSize',16);
end

%% Plot filters
%{
wML_full=full(wML);
selfmadeweight=zeros(1,926);
for i = 2:(nBases+1)
    selfmadeweight=selfmadeweight+wML_full(i)*bs.B(:,i-1)';
end
%}
wML_full = wML_full_both(1:6);
rateBias = wML_full(1);
w = zeros(length(bs.B(:,1)),1);
for i = 2:6
    w = w + wML_full(i)*bs.B(:,i-1);
end
if plotornot(2)==1
    figure; clf; hold all;
    %w=ws.hist.data;
    %w(1)=-Inf;
    %w(2)=-Inf;
    plot(w);
    %plot(selfmadeweight);
    %plot(est_filter);
    title('Input filter','FontSize',16);
    %legend('mine','original');
    xlabel(['t/ms'],'FontSize',13);
end

%% Plot filters
%{
wML_full=full(wML);
selfmadeweight=zeros(1,926);
for i = 2:(nBases+1)
    selfmadeweight=selfmadeweight+wML_full(i)*bs.B(:,i-1)';
end
%}
wML_full = wML_full_both(7:12);
rateBias = wML_full(1);
w = zeros(length(bs.B(:,1)),1);
for i = 2:6
    w = w + wML_full(i)*bs.B(:,i-1);
end
if plotornot(2)==1
    figure; clf; hold all;
    %w=ws.hist.data;
    %w(1)=-Inf;
    %w(2)=-Inf;
    plot(w);
    %plot(selfmadeweight);
    %plot(est_filter);
    title('History filter','FontSize',16);
    %legend('mine','original');
    xlabel(['t/ms'],'FontSize',13);
end


%% Plot KS test

if kswhich==0
    w=est_filter(2:end);
    filter_length=length(w);
end

%w(1)=-1e10;
%w(2)=-1e10;

fullsp=zeros(Ntrial,T);
for i=1:size(sp,1)
    for j=1:size(sp,2)
        if sp(i,j)~=0
            fullsp(i,sp(i,j))=1;
        end
    end
end

z=[];
rate=zeros(Ntrial,T);

for i=trainingtrials
    rate(i,1)=exp(rateBias);
    for t=2:T
        pasttime=min(filter_length,t-1);
        temp=w(1:pasttime)' * fullsp(i ,t - (1:pasttime))';
        rate(i,t) = exp( temp + rateBias );
    end
    z=[z,sum(rate(i,1:sp(i,1)))];
    for j=2:totalspikes(i)
        z=[z,sum(rate( i , (sp(i,j-1)+1) : sp(i,j) ))];
    end
	
end

if plotornot(3)==1
    figure
    [eCDF, zvals] = ecdf(z);
    mCDF = 1-exp(-zvals); 
    plot(mCDF,eCDF) 
    hold on 
    plot([0 1], [0 1]+1.36/sqrt(length(z)),'k') 
    plot([0 1], [0 1]-1.36/sqrt(length(z)),'k') 
    hold off 
    xlabel('Model CDF','FontSize',13) 
    ylabel('Empirical CDF','FontSize',13)
    title('Goodness-of-fit test','FontSize',16);
    axis([0 1 0 1]);
end

%% calulate P
%{
D=max(abs(mCDF-eCDF));
C=sqrt(length(raster_x))*D;
alpha=2*exp(-2*C^2);
display(alpha);
%}
test_cdf = makedist('exp', 'mu', 1);
[h, p] = kstest(z, 'CDF', test_cdf)

%% Plot simulated spike trains
%w(1)=-1e10;
%w(2)=-1e10;
if plotornot(4)==1
    simulationtimes=10;
    T=SimulationT;
    nHistBins = numel(w);
    y = zeros(simulationtimes , T + nHistBins);
    for NSimTrial=1:simulationtimes
        for t = (nHistBins+1):(T+nHistBins)
            yy = poissrnd(exp(w' * (y(NSimTrial,t - (1:nHistBins)))' + rateBias));
            if yy ~= 0
                y(NSimTrial,t) = 1;
            end
        end
    end

    y = y(:,nHistBins+1:end);
    st = find(y);

    figure
    plotraster(y,1:T,'Simulated Result');
end
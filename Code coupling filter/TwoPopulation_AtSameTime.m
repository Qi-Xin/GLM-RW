clearvars;
rng(1);

global hbasis Bspline nBspline nhbasis nPop T nTrial y

addpath(genpath('D:/code'))
nPop = 2;     % number of neuron population
nNeu = 1e2;    % number of neuons in a population
rec_nNeu = 1e0;      % number of neurons recorded in each population
T = 1e3;
nTrial = 3e2;
stopValue = 1e-3;
couplingStrength = 1/nNeu/3e1; % maximum of coupling filter
jitter = 0; % 0 for no jitter , 1 for standard amount of jitter
baselinefr = 3e-2;
highestfr = 1e-1;

maxIter = 10;
learningRate = 1; % change line 173 for slowly updating firing rate

dt = 1;
totT = nTrial*T;
fr = cell(1,nPop);
all_y = cell(1,nPop);
y = cell(1,nPop);

cp_spmodel = cell(nPop,nPop);
inhomoBias_spmodel = cell(1,nPop);
fr_spmodel = cell(1,nPop);

cp_frmodel = cell(nPop,nPop);
inhomoBias_frmodel = cell(1,nPop);
fr_frmodel = cell(1,nPop);
fr_frmodelNew = cell(1,nPop);

cp_true = cell(nPop,nPop);
%% Simulation
% Set fr1 and fr2
%baselinefr = 2e-2;

fr{1} = zeros(1,nTrial*T);
fr{2} = zeros(1,nTrial*T);
for i = 1:nTrial
    fr{1}( T*(i-1)+1 : T*i ) = get_signal(4,highestfr-baselinefr,T,jitter)+baselinefr;
%     fr{1}( T*(i-1)+1 : T*i ) = 3*baselinefr;
end

% get ground true coupling filter 1 and 2 (both are from a part of gamma)
temp = get_signal(3,couplingStrength,75,0);
cp1 = zeros(1,3e2);
cp1(1:length(temp)) = temp;

temp = get_signal(3,couplingStrength,3e2,0);
cp2 = zeros(1,3e2);
cp2(1:length(temp)) = temp;

cp_true{2,1} = cp1;
%cp_true{3,1} = cp2;

% get y1
all_y{1} = random('poisson',repmat(fr{1},nNeu,1));
%all_y{1} = random('binomial',1,repmat(fr{1},nNeu,1));

% get fr2
fr{2} = sameconv_Cutoff(sum(all_y{1})',cp1',T)'+baselinefr;
all_y{2} = random('poisson',repmat(fr{2},nNeu,1));
%all_y{2} = random('binomial',1,repmat(fr{2},nNeu,1));

% for i=1:nPop
%     all_y{i}(all_y{i}>1) = 1;
% end

% get observed spike trains
for i = 1:nPop
    y{i} = sum( all_y{i}(1:rec_nNeu,:),1 )';
end

%% spike train GLM

% make basis for post-spike kernel
ihbasprs.ncols = 7;  % number of basis vectors for post-spike kernel
hPeaksMax = 120;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 1*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms
[ht,hbas,hbasis] = makeBasis_PostSpike(ihbasprs,dt);
nhbasis = size(hbasis,2); % number of basis functions for h
% hbasis(2:end,1) = 0;

% make B spline basis for inhomogenerous underlying firing rate
nBspline = 3;
Bspline = makeBasis_spline(nBspline,T);
Bspline = repmat(Bspline,nTrial,1);

for i = 1:nPop
    yconvhi_all = [];
    for j = 1:nPop
        if i~=j
            yconvhi = zeros(size(y{j},1),nhbasis);
            for hnum = 1:nhbasis
                yconvhi(:,hnum) = sameconv(y{j},hbasis(:,hnum));
            end
            yconvhi_all = [yconvhi_all,yconvhi];
        end
    end
    [prs,dev,stats] = glmfit([Bspline,yconvhi_all],y{i},'poisson');
    se = stats.se;
    prs = [prs-se,prs,prs+se];
    dc = prs(1,:);
    B = Bspline*prs(2:nBspline+1,:);
    inhomoBias_spmodel{i} = B+dc;
    nn = 0;
    for j = 1:nPop
        if i~=j
            cp_spmodel{i,j} = hbasis*prs( (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1) , : );
            nn = nn+1;
        end
    end
    fr_spmodel{i} = exp( [Bspline,yconvhi_all]*prs(2:end,2)+prs(1,2) );
end

nlogL = 0;
for i = 1:nPop
    nlogL = nlogL + sum( fr_spmodel{i}-y{i}.*log(fr_spmodel{i}) ) ;
end
nlogL_spmodel = nlogL

%% firing rate GLM

Bspline = makeBasis_spline(nBspline,T);
fun = @NLogLikelihood;
x0 = 1e-2*ones(1,2*(1+nhbasis+nBspline));
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','MaxIterations',100, ...
    'OptimalityTolerance',1e-10);
[x,fval] = fminunc(fun,x0,options)

[fr_frmodel,cp_frmodel,inhomoBias_frmodel] = getFR(x);

%%
figure
for i = 1:nPop
    for j = 1:nPop
        if i~=j
            subplot(4,3,3*i-3+j)
            hold on
            plot(cp_spmodel{i,j},'-b','LineWidth',1);
            %plot(cp_frmodel{i,j},'-r','LineWidth',1);
            if isempty(cp_true{i,j})
                plot(zeros(1,200),'-k','LineWidth',1);
            else
                plot(cp_true{i,j},'-k','LineWidth',1);
            end
        end
    end
end
for i = 1:nPop
    subplot(4,3,9+i)
    hold on
    plot(inhomoBias_spmodel{i}(1:T,:),'-b','LineWidth',1);
    %plot(inhomoBias_frmodel{i}(1:T),'-r','LineWidth',1);
end

%%
figure
for i = 1:nPop
    for j = 1:nPop
        if i~=j
            subplot(4,3,3*i-3+j)
            hold on
            %plot(cp_spmodel{i,j},'-b','LineWidth',1);
            plot(cp_frmodel{i,j},'-r','LineWidth',1);
            if isempty(cp_true{i,j})
                plot(zeros(1,200),'-k','LineWidth',1);
            else
                plot(cp_true{i,j},'-k','LineWidth',1);
            end
        end
    end
end
for i = 1:nPop
    subplot(4,3,9+i)
    hold on
    %plot(inhomoBias_spmodel{i}(1:T,:),'-b','LineWidth',1);
    plot(inhomoBias_frmodel{i}(1:T),'-r','LineWidth',1);
end
for i = 1:nPop
    subplot(4,3,9+i)
    hold on
    %plot(inhomoBias_spmodel{i}(1:T,:),'-b','LineWidth',1);
    plot(inhomoBias_frmodel{i}(1:T),'-r','LineWidth',1);
end
%%
figure
subplot(2,1,1)
plot(prs_rcd{1}');
legend('\beta_1^{(i)}','\beta_2^{(i)}','\beta_3^{(i)}','\beta_4^{(i)}','\beta_5^{(i)}');
title('Coefficients for the Coupling filter from Neuron 1 to Neuron 2');
xlabel('Iteration i');
subplot(2,1,2)
plot(prs_rcd{2}');
legend('\beta_1^{(i)}','\beta_2^{(i)}','\beta_3^{(i)}','\beta_4^{(i)}','\beta_5^{(i)}');
xlabel('Iteration i');
title('Coefficients for the Coupling filter from Neuron 2 to Neuron 1');

%% function defination
function [fr_frmodel,cp_frmodel,inhomoBias_frmodel] = getFR(x)
    global hbasis Bspline nBspline nhbasis nPop T nTrial y
    nprs = (1+nBspline+nhbasis*(nPop-1));
    cp_frmodel = cell(nPop,nPop);
    inhomoBias_frmodel = cell(1,nPop);
    fr_frmodel = cell(1,nPop);
    for i = 1:nPop
        prs{i} = x(1+(i-1)*nprs:i*nprs)';
        dc{i} = prs{i}(1);
        B{i} = Bspline*prs{i}(2:nBspline+1);
        inhomoBias_frmodel{i} = B{i} + dc{i};
        fr_frmodel{i} = inhomoBias_frmodel{i};

        nn = 0;
        for j = 1:nPop
            if i~=j
                cp_frmodel{i,j} = hbasis*prs{i}( (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1) , : );
                nn = nn+1;
            end
        end
    end
	
    
    for t = 2:T
        for i = 1:nPop
            for j = 1:nPop
                if i~=j
                    FRNow = sameconv( fr_frmodel{j}(1:t-1) , cp_frmodel{i,j});
                    fr_frmodel{i}(t) = fr_frmodel{i}(t) + FRNow(end);
                end
            end
        end
    end
    
end


function output_args = NLogLikelihood(x)
    global hbasis Bspline nBspline nhbasis nPop T nTrial y
    [fr_frmodel,cp_frmodel,inhomoBias_frmodel] = getFR(x);
    for i = 1:nPop
        fr_frmodel{i} = repmat(fr_frmodel{i},nTrial,1);
    end
    nlogL = 0;
    for i = 1:nPop
        nlogL = nlogL + sum( fr_frmodel{i}-y{i}.*log(fr_frmodel{i}) ) ;
    end
    output_args = nlogL;
    
end
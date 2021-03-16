clearvars;
rng(1);

addpath(genpath('D:/code'))
nPop = 2;     % number of neuron population
nNeu = 1e2;    % number of neuons in a population
rec_nNeu = 1e0;      % number of neurons recorded in each population
T = 1e3;
nTrial = 5e2;
stopValue = 1e-3;
couplingStrength = 1/nNeu/1e2; % maximum of coupling filter
jitter = 0;

dt = 1;
totT = nTrial*T;
fr = cell(1,nPop);
all_y = cell(1,nPop);
y = cell(1,nPop);
y_smooth = cell(1,nPop);

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
baselinefr = 1e-2;
highestfr = 7e-2;

fr{1} = zeros(1,nTrial*T);
fr{2} = zeros(1,nTrial*T);
for i = 1:nTrial
    fr{1}( T*(i-1)+1 : T*i ) = get_signal(3,highestfr-baselinefr,T,jitter)+baselinefr;
end

% get ground true coupling filter 1 and 2 (both are from a part of gamma)
temp = get_signal(3,couplingStrength,75,0);
temp = temp(5:end);
cp1 = zeros(1,3e2);
cp1(1:length(temp)) = temp;

temp = get_signal(3,couplingStrength,3e2,0);
temp = temp(30:end);
cp2 = zeros(1,3e2);
cp2(1:length(temp)) = temp;

cp_true{2,1} = cp1;
%cp_true{3,1} = cp2;

% get y1
all_y{1} = random('poisson',repmat(fr{1},nNeu,1));
%all_y{1} = random('binomial',1,repmat(fr{1},nNeu,1));

% get fr2
fr{2} = sameconv(sum(all_y{1})',cp1')'+baselinefr;
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
ihbasprs.ncols = 2;  % number of basis vectors for post-spike kernel
hPeaksMax = 40;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 1e3*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms
[ht,hbas,hbasis] = makeBasis_PostSpike(ihbasprs,dt);
nhbasis = size(hbasis,2); % number of basis functions for h
% hbasis(2:end,1) = 0;

% make B spline basis for inhomogenerous underlying firing rate
nBspline = 2;
Bspline = makeBasis_spline(nBspline,T);
Bspline = repmat(Bspline,nTrial,1);

smoothSigmaList = ceil(logspace(0,3,30));
nlogL_spmodel= [];
for smoothSigma = smoothSigmaList
    for i = 1:nPop
        xx = [-3*smoothSigma:1:3*smoothSigma];
        yy = normpdf(xx,0,smoothSigma);
        y_smooth{i} = conv(y{i},yy);
        y_smooth{i} = y_smooth{i}(3*smoothSigma+1:end-3*smoothSigma);
    end
    for i = 1:nPop
        yconvhi_all = [];
        for j = 1:nPop
            if i~=j
                yconvhi = zeros(size(y{j},1),nhbasis);
                for hnum = 1:nhbasis
                    yconvhi(:,hnum) = sameconv(y_smooth{j},hbasis(:,hnum));
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
    nlogL_spmodel = [nlogL_spmodel,nlogL];
end
%%
figure;
plot(smoothSigmaList,nlogL_spmodel);

%%
smoothSigma = smoothSigmaList(find(nlogL_spmodel == min(nlogL_spmodel)));
figure
for i = 1:nPop
    xx = [-3*smoothSigma:1:3*smoothSigma];
    yy = normpdf(xx,0,smoothSigma);
    y_smooth{i} = conv(y{i},yy);
    y_smooth{i} = y_smooth{i}(3*smoothSigma+1:end-3*smoothSigma);
end
for i = 1:nPop
    yconvhi_all = [];
    for j = 1:nPop
        if i~=j
            yconvhi = zeros(size(y{j},1),nhbasis);
            for hnum = 1:nhbasis
                yconvhi(:,hnum) = sameconv(y_smooth{j},hbasis(:,hnum));
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
nlogL_spmodelBest = nlogL;

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
title('sp model');

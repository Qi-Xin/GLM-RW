clearvars;

addpath(genpath('D:/code'))
nPop = 3;     % number of neuron population
nNeu = 1e2;    % number of neuons in a population
rec_nNeu = 1;      % number of neurons recorded in each population
T = 5e3;
nTrial = 1e2;
stopValue = 1e-2;
maxIter = 1e3;
couplingStrength = 1/nNeu/3e2; % maximum of coupling filter

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
baselinefr = 2e-2;
highestfr = 7e-2;
fr{1} = get_signal(3,highestfr-baselinefr,T)+baselinefr;
fr{2} = get_signal(4,highestfr-baselinefr,T)+baselinefr;

% get ground true coupling filter 1 and 2 (both are from a part of gamma)
temp = get_signal(3,couplingStrength,2.5e2);
temp = temp(15:end);
cp1 = zeros(1,3e2);
cp1(1:length(temp)) = temp;

temp = get_signal(3,couplingStrength,3e2);
temp = temp(30:end);
cp2 = zeros(1,3e2);
cp2(1:length(temp)) = temp;

cp_true{2,1} = cp1;
cp_true{3,1} = cp2;

% get y1 and y2
all_y{1} = random('poisson',repmat(fr{1},nNeu,nTrial));
all_y{2} = random('poisson',repmat(fr{2},nNeu,nTrial));

% get fr3
fr{3} = sameconv(sum(all_y{1})',cp1')'+sameconv(sum(all_y{2})',cp2')'+baselinefr;
all_y{3} = random('poisson',repmat(fr{3},nNeu,1));

for i=1:nPop
    all_y{i}(all_y{i}>1) = 1;
end

% get observed spike trains
for i = 1:nPop
    y{i} = all_y{i}(1,:)';
end
%%
i = 1;
figure
subplot(2,1,1);
plotraster(reshape(y{i}(1:nTrial*T),[],nTrial)',1:T,'Simulated Result');
title('Spike train');
xlabel('ms');
ylabel('trial');
%% spike train GLM
% make basis for post-spike kernel
ihbasprs.ncols = 5;  % number of basis vectors for post-spike kernel
hPeaksMax = 75;
ihbasprs.hpeaks = [0 hPeaksMax];  % peak location for first and last vectors, in ms
ihbasprs.b = 0.2*hPeaksMax;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms
[ht,hbas,hbasis] = makeBasis_PostSpike(ihbasprs,dt);
nhbasis = size(hbasis,2); % number of basis functions for h

% make B spline basis for inhomogenerous underlying firing rate
nBspline = 5;
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
    
    if i==3
        [prs,dev,stats] = glmfit([yconvhi_all],y{i},'poisson');
        se = stats.se;
        prs = [prs-se,prs,prs+se];
        dc = prs(1,:);
        inhomoBias_spmodel{i} = dc;
        nn = 0;
        for j = 1:nPop
            if i~=j
                cp_spmodel{i,j} = hbasis*prs( (nhbasis*nn+2):(nhbasis*(nn+1)+1) , : );
                nn = nn+1;
            end
        end
        fr_spmodel{i} = exp( [yconvhi_all]*prs(2:end,:) );
    else
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
    

    
end


%% firing rate GLM

fr_frmodel = y;
nlogLBest = Inf;
noImproveIter = 0;
for iter = 1:maxIter
    iter
    for i = 1:nPop
        yconvhi_all = [];
        for j = 1:nPop
            if i~=j
                yconvhi = zeros(size(fr_frmodel{j},1),nhbasis);
                for hnum = 1:nhbasis
                    yconvhi(:,hnum) = sameconv(fr_frmodel{j},hbasis(:,hnum));
                end
                yconvhi_all = [yconvhi_all,yconvhi];
            end
        end
        [prs,dev,stats] = glmfit([Bspline,yconvhi_all],y{i},'poisson');
        se = stats.se;
        prs = [prs-se,prs,prs+se];
        fr_frmodelNew{i} = exp( [Bspline,yconvhi_all]*prs(2:end,2)+prs(1,2) );
        dc = prs(1,:);
        B = Bspline*prs(2:nBspline+1,:);
        inhomoBias_frmodel{i} = B+dc;
        nn = 0;
        for j = 1:nPop
            if i~=j
                cp_frmodel{i,j} = hbasis*prs( (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1) , : );
                nn = nn+1;
            end
        end
        
        diff = max(fr_frmodel{i}-fr_frmodelNew{i})
        if diff <= stopValue
            noImproveIter = noImproveIter+1;
            if noImproveIter >= 3
                break
            end
        else
            noImproveIter = 0;
        end
        
        fr_frmodel{i} = fr_frmodelNew{i};

        
    end
    
    if noImproveIter >= 3
        break
    end
    
%     nlogL = 0;
%     for i = 1:nPop
%         nlogL = nlogL + sum( fr_frmodel{i}-y{i}.*log(fr_frmodel{i}) ) ;
%     end
%     if isnan(nlogL)
%         nlogL = Inf;
%     end
%     if nlogL <= nlogLBest
%         cp_frmodelBest = cp_frmodel;
%         nlogLBest = nlogL;
%         noImproveIter = 0;
%     else
%         noImproveIter = noImproveIter+1;
%         if noImproveIter >= 3
%             break
%         end
%     end
%     display([nlogL,nlogLBest,noImproveIter]);

end


%%
figure
for i = 1:nPop
    for j = 1:nPop
        if i~=j
            subplot(5,3,3*i-3+j)
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
for i = 1:2
    subplot(5,3,9+i)
    hold on
    plot(inhomoBias_spmodel{i}(1:T,:),'-b','LineWidth',1);
    plot(inhomoBias_frmodel{i}(1:T,:),'-r','LineWidth',1);
end
for i = 1:nPop
    subplot(5,3,12+i)
    hold on
    plot(fr_spmodel{i}(1:T,:),'-b','LineWidth',1);
    plot(inhomoBias_frmodel{i}(1:T,:),'-r','LineWidth',1);
end
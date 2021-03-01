clearvars;
rng(1);

addpath(genpath('D:/code'))
nPop = 2;     % number of neuron population
nNeu = 1e2;    % number of neuons in a population
rec_nNeu = 1e0;      % number of neurons recorded in each population
T = 1e3;
nTrial = 1e0;
stopValue = 1e-3;
maxIter = 10;
couplingStrength = 1/nNeu/3e2; % maximum of coupling filter
jitter = 0;

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
baselinefr = 1e-2;
highestfr = 7e-2;

fr{1} = zeros(1,nTrial*T);
fr{2} = zeros(1,nTrial*T);
for i = 1:nTrial
    fr{1}( T*(i-1)+1 : T*i ) = 1e-2*ones(1,T);
end

% get ground true coupling filter 1 and 2 (both are from a part of gamma)
temp = get_signal(3,couplingStrength,75,0);
temp = temp(15:end);
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

% make B spline basis for inhomogenerous underlying firing rate
nBspline = 2;
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
    fr_spmodel{i} = exp( [Bspline,yconvhi_all]*prs(2:end,:) );
end


%% firing rate GLM
prs_rcd = cell(1,2);
nlogL_rcd = [];
maxIter = 50;
color = {'r','g','b'};
figure

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
        prs_rcd{i} = [prs_rcd{i},prs];
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
        
        if iter>47
            subplot(2,2,2*i)
            hold on
            temp = cp_frmodel{i,3-i}(:,2);
            plot(temp,color{mod(iter,3)+1});
            text(0,temp(1),[num2str(iter),num2str(i)]);
            
            subplot(2,2,2*i-1)
            hold on
            
            yyaxis left;
            temp = fr_frmodel{i}(:,1);
            plot(temp,color{mod(iter,3)+1});
            text(0,temp(1),[num2str(iter),num2str(i)]);
            yyaxis right;
            plotraster(y{i}',1:T,'Simulated Result');
            ylim([0 4]);
        end

        
    end
    
    nlogL = 0;
    for i = 1:nPop
        nlogL = nlogL + sum( fr_frmodel{i}-y{i}.*log(fr_frmodel{i}) ) ;
    end
    nlogL_rcd = [nlogL_rcd,nlogL];
    
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
figure;
hold on;
yyaxis left;
plot(fr_frmodel{1}(:,1))

yyaxis right;
plotraster(y{1}',1:T,'Simulated Result');
ylim([0 10]);

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
subplot(3,1,1)
plot(prs_rcd{1}');
subplot(3,1,2)
plot(prs_rcd{2}');
subplot(3,1,3)
plot(nlogL_rcd');
clearvars;
rng(1);

addpath(genpath('D:/Github/GLM-RW'))
nPop = 2;     % number of neuron population
nNeu = 1e2;    % number of neuons in a population
rec_nNeu = 1e0;      % number of neurons recorded in each population
T = 1e3;
nTrial = 3e2;
stopValue = 1e-3;
couplingStrength = 1/nNeu/3e1; % maximum of coupling filter
jitter = 1; % 0 for no jitter , 1 for standard amount of jitter
baselinefr = 3e-2;
highestfr = 1e-1;

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
                    yconvhi(:,hnum) = sameconv_Cutoff(y_smooth{j},hbasis(:,hnum),T);
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
smoothSigma = min( smoothSigmaList(find(nlogL_spmodel == min(nlogL_spmodel))) , 200) ;
%%
figure(1);
plot(smoothSigmaList,nlogL_spmodel);
%% Best smooth para

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
                yconvhi(:,hnum) = sameconv_Cutoff(y_smooth{j},hbasis(:,hnum),T);
            end
            yconvhi_all = [yconvhi_all,yconvhi];
        end
    end
    [prs,dev,stats] = glmfit([Bspline,yconvhi_all],y{i},'poisson');
    se = stats.se;
    covb = stats.covb;
    prs = [prs-se,prs,prs+se];
    dc = prs(1,:);
    B = Bspline*prs(2:nBspline+1,:);
    inhomoBias_spmodel{i} = B+dc;
    nn = 0;
    for j = 1:nPop
        if i~=j
            cp_spmodel{i,j} = hbasis*prs( (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1) , : );
            cov_cp = covb((nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1),...
                (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1));
            var_of_cp = hbasis*cov_cp*hbasis';
            cp_spmodel{i,j}(:,1) = cp_spmodel{i,j}(:,2) - sqrt(diag(var_of_cp));
            cp_spmodel{i,j}(:,3) = cp_spmodel{i,j}(:,2) + sqrt(diag(var_of_cp));
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


figure(2)
nplot = 0;
for i = 1:nPop
    for j = 1:nPop
        if i~=j
            nplot = nplot+1;
            subplot(3,2,nplot)
            k = cp_spmodel{i,j}(:,2);
            kL = cp_spmodel{i,j}(:,1);
            kU = cp_spmodel{i,j}(:,3);
            fill([(1:length(k)) fliplr(1:length(k))],[kL' fliplr(kU')], ...
                'b','facealpha',0.2,'edgealpha',0,'HandleVisibility','off');
            hold on
            plot(k,'-b','LineWidth',1);
            if isempty(cp_true{i,j})
                plot(zeros(1,200),'-r','LineWidth',1.5);
            else
                plot(cp_true{i,j}/max(cp_true{i,j})*max(cp_spmodel{i,j}(:,2)),'-r','LineWidth',1.5);
            end
            legend('Fitted ','Ground True');
            xlim([0 200])
            title('Coupling Filter from Neuron A to Neuron B');
        end
    end
end
sample_estimated_cp_sm = cp_spmodel{2,1};
sample_estimated_fr_sm = fr_spmodel{2};
sample_estimated_cppart_sm = yconvhi_all*prs( (nBspline+2):(nBspline+nhbasis+1) , 2 );
%% Raw spike train
for i = 1:nPop
    yconvhi_all = [];
    for j = 1:nPop
        if i~=j
            yconvhi = zeros(size(y{j},1),nhbasis);
            for hnum = 1:nhbasis
                yconvhi(:,hnum) = sameconv_Cutoff(y{j},hbasis(:,hnum),T);
            end
            yconvhi_all = [yconvhi_all,yconvhi];
        end
    end
    [prs,dev,stats] = glmfit([Bspline,yconvhi_all],y{i},'poisson');
    se = stats.se;
    covb = stats.covb;
    prs = [prs-se,prs,prs+se];
    dc = prs(1,:);
    B = Bspline*prs(2:nBspline+1,:);
    inhomoBias_spmodel{i} = B+dc;
    nn = 0;
    for j = 1:nPop
        if i~=j
            cp_spmodel{i,j} = hbasis*prs( (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1) , : );
            cov_cp = covb((nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1),...
                (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1));
            var_of_cp = hbasis*cov_cp*hbasis';
            cp_spmodel{i,j}(:,1) = cp_spmodel{i,j}(:,2) - sqrt(diag(var_of_cp));
            cp_spmodel{i,j}(:,3) = cp_spmodel{i,j}(:,2) + sqrt(diag(var_of_cp));
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

figure(2)
for i = 1:nPop
    for j = 1:nPop
        if i~=j
            nplot = nplot+1;
            subplot(3,2,nplot)
            k = cp_spmodel{i,j}(:,2);
            kL = cp_spmodel{i,j}(:,1);
            kU = cp_spmodel{i,j}(:,3);
            fill([(1:length(k)) fliplr(1:length(k))],[kL' fliplr(kU')], ...
                'b','facealpha',0.2,'edgealpha',0,'HandleVisibility','off');
            hold on
            plot(k,'-b','LineWidth',1);
            if isempty(cp_true{i,j})
                plot(zeros(1,200),'-r','LineWidth',1.5);
            else
                plot(cp_true{i,j}/max(cp_true{i,j})*max(cp_spmodel{i,j}(:,2)),'-r','LineWidth',1.5);
            end
            legend('Fitted ','Ground True');
            xlim([0 200])
            title('Coupling Filter from Neuron A to Neuron B');
        end
    end
end
sample_estimated_cp_raw = cp_spmodel{2,1};
sample_estimated_fr_raw = fr_spmodel{2};
sample_estimated_cppart_raw = yconvhi_all*prs( (nBspline+2):(nBspline+nhbasis+1) , 2 );
%% Pooling spike train
for i = 1:nPop
    yconvhi_all = [];
    for j = 1:nPop
        if i~=j
            yconvhi = zeros(size(y{j},1),nhbasis);
            for hnum = 1:nhbasis
                pooling = sum(reshape(y{j},T,[])')';
                pooling = repmat(pooling,length(y{j})/length(pooling),1);
                yconvhi(:,hnum) = sameconv_Cutoff(pooling,hbasis(:,hnum),T);
            end
            yconvhi_all = [yconvhi_all,yconvhi];
        end
    end
    [prs,dev,stats] = glmfit([Bspline,yconvhi_all],y{i},'poisson');
    se = stats.se;
    covb = stats.covb;
    prs = [prs-se,prs,prs+se];
    dc = prs(1,:);
    B = Bspline*prs(2:nBspline+1,:);
    inhomoBias_spmodel{i} = B+dc;
    nn = 0;
    for j = 1:nPop
        if i~=j
            cp_spmodel{i,j} = hbasis*prs( (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1) , : );
            cov_cp = covb((nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1),...
                (nBspline+nhbasis*nn+2):(nBspline+nhbasis*(nn+1)+1));
            var_of_cp = hbasis*cov_cp*hbasis';
            cp_spmodel{i,j}(:,1) = cp_spmodel{i,j}(:,2) - sqrt(diag(var_of_cp));
            cp_spmodel{i,j}(:,3) = cp_spmodel{i,j}(:,2) + sqrt(diag(var_of_cp));
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

figure(2)
for i = 1:nPop
    for j = 1:nPop
        if i~=j
            nplot = nplot+1;
            subplot(3,2,nplot)
            k = cp_spmodel{i,j}(:,2);
            kL = cp_spmodel{i,j}(:,1);
            kU = cp_spmodel{i,j}(:,3);
            fill([(1:length(k)) fliplr(1:length(k))],[kL' fliplr(kU')], ...
                'b','facealpha',0.2,'edgealpha',0,'HandleVisibility','off');
            hold on
            plot(k,'-b','LineWidth',1);
            if isempty(cp_true{i,j})
                plot(zeros(1,200),'-r','LineWidth',1.5);
            else
                plot(cp_true{i,j}/max(cp_true{i,j})*max(cp_spmodel{i,j}(:,2)),'-r','LineWidth',1.5);
            end
            legend('Fitted ','Ground True');
            xlim([0 200])
            title('Coupling Filter from Neuron A to Neuron B');
        end
    end
end

    
sample_estimated_cp_pooling = cp_spmodel{2,1};
sample_estimated_fr_pooling = fr_spmodel{2};
sample_estimated_cppart_pooling = yconvhi_all*prs( (nBspline+2):(nBspline+nhbasis+1) , 2 );

%%
theTrialToLookAt = 5;
theTimeToLookAt = (theTrialToLookAt*T-T+301):(theTrialToLookAt*T-T+700);
figure(4)
subplot(4,3,1)
GroundTrueSpTrain = all_y{1}(:,theTimeToLookAt);
plotraster(GroundTrueSpTrain,1:400,'Simulated Result');
title('Ground true input spike train');
subplot(4,3,2)
plot(cp_true{2,1},'-k','LineWidth',1.5);
title('Ground true coupling filter');
subplot(4,3,3)
plot(theTimeToLookAt,fr{2}(theTimeToLookAt),'-k','LineWidth',1.5);
title('Ground true firing rate');

subplot(4,3,4)
SpTrain = y{1}(theTimeToLookAt);
plotraster(SpTrain',1:400,'Simulated Result');
title('Raw spike train recorded');
subplot(4,3,5)
k = sample_estimated_cp_raw(:,2);
kL = sample_estimated_cp_raw(:,1);
kU = sample_estimated_cp_raw(:,3);
fill([(1:length(k)) fliplr(1:length(k))],[kL' fliplr(kU')], ...
    'b','facealpha',0.2,'edgealpha',0,'HandleVisibility','off');
hold on
plot(k,'-b','LineWidth',1);
title('Raw spike train estimation of coupling filter ');
subplot(4,3,6)
plot(theTimeToLookAt,sample_estimated_fr_raw(theTimeToLookAt),'-b','LineWidth',1.5);
title('Raw spike train estimation of firing rate');

subplot(4,3,7)
plot(theTimeToLookAt,y_smooth{1}(theTimeToLookAt),'-k','LineWidth',1.5);
title('Smoothed spike train');
subplot(4,3,8)
k = sample_estimated_cp_sm(:,2);
kL = sample_estimated_cp_sm(:,1);
kU = sample_estimated_cp_sm(:,3);
fill([(1:length(k)) fliplr(1:length(k))],[kL' fliplr(kU')], ...
    'r','facealpha',0.2,'edgealpha',0,'HandleVisibility','off');
hold on
plot(k,'-r','LineWidth',1);
title('Smoothed spike train estimation of coupling filter ');
subplot(4,3,9)
plot(theTimeToLookAt,sample_estimated_fr_sm(theTimeToLookAt),'-r','LineWidth',1.5);
title('Smoothed spike train estimation of firing rate');

subplot(4,3,10)
plot(theTimeToLookAt,pooling(theTimeToLookAt),'-k','LineWidth',1.5);
title('Pooled spike train');
subplot(4,3,11)
k = sample_estimated_cp_pooling(:,2);
kL = sample_estimated_cp_pooling(:,1);
kU = sample_estimated_cp_pooling(:,3);
fill([(1:length(k)) fliplr(1:length(k))],[kL' fliplr(kU')], ...
    'g','facealpha',0.2,'edgealpha',0,'HandleVisibility','off');
hold on
plot(k,'-g','LineWidth',1);
title('Pooled spike train estimation of coupling filter ');
subplot(4,3,12)
plot(theTimeToLookAt,sample_estimated_fr_pooling(theTimeToLookAt),'-g','LineWidth',1.5);
title('Pooled spike train estimation of firing rate');

%%
theTrialToLookAt = 5;
theTimeToLookAt = (theTrialToLookAt*T-T+301):(theTrialToLookAt*T-T+700);
figure(3)
subplot(3,1,1)
hold on
plot(theTimeToLookAt,fr{1}(theTimeToLookAt));
plot(theTimeToLookAt,y_smooth{1}(theTimeToLookAt));
plot(theTimeToLookAt,y{1}(theTimeToLookAt)/30);
plot(theTimeToLookAt,pooling(theTimeToLookAt)/600);
subplot(3,1,2)
hold on
plot(theTimeToLookAt,fr{2}(theTimeToLookAt));
plot(theTimeToLookAt,sample_estimated_fr_sm(theTimeToLookAt));
plot(theTimeToLookAt,sample_estimated_fr_raw(theTimeToLookAt));
plot(theTimeToLookAt,sample_estimated_fr_pooling(theTimeToLookAt));
subplot(3,1,3)
hold on
plot(theTimeToLookAt,fr{2}(theTimeToLookAt)-baselinefr);
plot(theTimeToLookAt,sample_estimated_cppart_sm(theTimeToLookAt));
plot(theTimeToLookAt,sample_estimated_cppart_raw(theTimeToLookAt));
plot(theTimeToLookAt,sample_estimated_cppart_pooling(theTimeToLookAt));


legend('Ground true','sm','raw','pooling');
function [k h dc prs kbasis hbasis stats] = fit_FNL(x,y,dt,nkt,kbasprs,ihbasprs,fit_k,plotFlag,NumSP,ISI)
% [k h dc prs kbasis hbasis] = fit_glm(x,y,dt,nkt,kbasprs,ihbasprs,prs,softRect,plotFlag,maxIter,tolFun,L2pen)
%
%  This code fits a Poisson GLM to given data, using basis vectors to
%  characterize the stimulus and post-spike filters.
%
%  The inputs are:
%   x: stimulus
%   y: spiking data, vector of 0s and 1s
%   dt: time step of x and y in ms
%   nkt: number of ms in stimulus filter
%   kbasprs: structure containing parameters of stimulus filter basis vectors
%       kbasprs.neye: number of "identity" basis vectors near time of spike
%       kbasprs.ncos: number of raised-cosine vectors to use
%       kbasprs.kpeaks: position of first and last bump (relative to identity bumps)
%       kbasprs.b: how nonlinear to make spacings (larger -> more linear)
%   ihbasprs: structure containing parameters of post-spike filter basis vectors
%       ihbasprs.ncols: number of basis vectors for post-spike kernel
%       ihbasprs.hpeaks: peak location for first and last vectors
%       ihbasprs.b: how nonlinear to make spacings (larger -> more linear)
%       ihbasprs.absref: absolute refractory period, in ms
%   prs: vector to initialize fit parameters
%   softRect: 0 uses exponential nonlinearity; 1 uses soft-rectifying nonlinearity
%   plotFlag: 0 or 1, plot simulated data
%   maxIter: maximum number of iterations for fitting
%   tolFun: function tolerance for fitting
%   L2pen: size of L2 penalty on coefficients in prs (defaults to 0)
%
%  The outputs are:
%   k: stimulus filter
%   h: post-spike filter
%   dc: DC offset
%   prs: full set of coefficients for basis vectors, [k_coeffs h_coeffs dc]
%   kbasis: basis vectors for stimulus filter
%   hbasis: basis vectors for post-spike filters
%
%  This code requires the function fminunc in the MATLAB Optimization
%  Toolbox, as well as the following functions:
%   makeBasis_StimKernel
%   makeBasis_PostSpike
%   normalizecols
%   sameconv
%   negloglike_glm_basis (or negloglike_glm_basis_softrect, if using soft-rectified nonlinearity)
%   logexp1 (if using soft-rectified nonlinearity)

%% set defaults

if ~exist('nkt','var') || isempty(nkt)
    nkt = 100;
end

if ~exist('kbasprs','var') || isempty(kbasprs)
    %%% basis functions for stimulus filter
    kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
    kbasprs.ncos = 3; % number of raised-cosine vectors to use
    kbasprs.kpeaks = [1 round(nkt/2)];  % position of first and last bump (relative to identity bumps)
    kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
end

if ~exist('ihbasprs','var') || isempty(ihbasprs)
    %%% basis functions for post-spike kernel
    ihbasprs.ncols = 2;  % number of basis vectors for post-spike kernel
    ihbasprs.hpeaks = [1 100];  % peak location for first and last vectors,in ms
    ihbasprs.b = 10;  % how nonlinear to make spacings (larger -> more linear)
    ihbasprs.absref = 0; % absolute refractory period, in ms
end

if ~exist('softRect','var') || isempty(softRect)
    softRect = 0;
end

if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = 0;
end

if ~exist('maxIter','var') || isempty(maxIter)
    maxIter = 100;
end

if ~exist('tolFun','var') || isempty(tolFun)
    tolFun = 1e-8;
end

if ~exist('L2pen','var') || isempty(L2pen)
    L2pen = 0;  % penalty on L2 norm
end

refreshRate = 1000/dt; % stimulus in ms, sampled at dt


%% create basis functions and initialize parameters

kbasisTemp = makeBasis_StimKernel(kbasprs,nkt);
nkb = size(kbasisTemp,2);
lenkb = size(kbasisTemp,1);
kbasis = zeros(lenkb/dt,nkb);
for bNum = 1:nkb
    kbasis(:,bNum) = interp1([1:lenkb]',kbasisTemp(:,bNum),linspace(1,lenkb,lenkb/dt)');
end

[ht,hbas,hbasis] = makeBasis_PostSpike(ihbasprs,dt);
hbasis = [zeros(1,ihbasprs.ncols); hbasis]; % enforce causality: post-spike filter only affects future time points

nkbasis = size(kbasis,2); % number of basis functions for k
nhbasis = size(hbasis,2); % number of basis functions for h

if ~exist('prs','var') || isempty(prs)
    prs = zeros(nkbasis+nhbasis+1,1); % initialize parameters
end

%%

xconvki = zeros(size(y,1),nkbasis);
yconvhi = zeros(size(y,1),NumSP*nhbasis);

for knum = 1:nkbasis
    xconvki(:,knum) = sameconvSti(x,flipud(kbasis(:,knum)));
end

for i = 1:NumSP
    for hnum = 1:nhbasis
        yconvhi(:,hnum+(i-1)*nhbasis) = sameconv_FNL(y,hbasis(:,hnum),i);
    end
end

%% optimization with glmfit
if fit_k == 1
    [prs,dev,stats] = glmfit([xconvki,yconvhi],y,'poisson');
else
    [prs,dev,stats] = glmfit([yconvhi],y,'poisson');
end
se = stats.se;
prsL = prs-se;
prsU= prs+se;
h = cell(1,NumSP);
hL = cell(1,NumSP);
hU = cell(1,NumSP);
%% calculate filters from basis fcns/weights
if fit_k == 1
    k = kbasis*prs(2:nkbasis+1); % k basis functions weighted by given parameters
    kL = kbasis*prsL(2:nkbasis+1);
    kU = kbasis*prsU(2:nkbasis+1);
    for i = 1:NumSP
        h{i} = hbasis*prs(nkbasis+(i-1)*nhbasis+2:nkbasis+i*nhbasis+1); % k basis functions weighted by given parameters
        hL{i} = hbasis*prsL(nkbasis+(i-1)*nhbasis+2:nkbasis+i*nhbasis+1);
        hU{i} = hbasis*prsU(nkbasis+(i-1)*nhbasis+2:nkbasis+i*nhbasis+1);
        h{i} = h{i}(1:end);
        hL{i} = hL{i}(1:end);
        hU{i} = hU{i}(1:end);
    end
else
    k = kbasis*zeros(nkbasis,1);
    kL = kbasis*zeros(nkbasis,1);
    kU = kbasis*zeros(nkbasis,1);
    for i = 1:NumSP
        h{i} = hbasis*prs((i-1)*nhbasis+2:i*nhbasis+1); % k basis functions weighted by given parameters
        hL{i} = hbasis*prsL((i-1)*nhbasis+2:i*nhbasis+1);
        hU{i} = hbasis*prsU((i-1)*nhbasis+2:i*nhbasis+1);
        h{i} = h{i}(1:end);
        hL{i} = hL{i}(1:end);
        hU{i} = hU{i}(1:end);
    end
end
dc = prs(1); % dc current (accounts for mean spike rate)
dcLU = [prsL(1),prsU(1)];

%% plot results
if plotFlag
    figure;

    subplot(3,1,1);hold on
    plot((1:length(k)),k);
    if fit_k == 1
        fill([(1:length(k)) fliplr(1:length(k))],[kL' fliplr(kU')],'b','facealpha',0.2,'edgealpha',0);
    end
    xlim([0 length(k) ]);
    set(gca,'xtick',0:round(length(k)/2):length(k),'xticklabel',round(fliplr(-dt*(0:round(length(k)/2):length(k)))))
    xlabel('time (ms)')
    title('stimulus filter')
    grid on
    
    subplot(3,1,3);hold on
    for i = 1:NumSP
        plot((1:length(h{i})),h{i});
        fill([(1:length(h{i})) fliplr(1:length(h{i}))],[hL{i}' fliplr(hU{i}')],'b','facealpha',0.2,'edgealpha',0,'HandleVisibility','off');
    end
    xlim([1 length(h{i}) ]);
    legend('1','2','3');
    %set(gca,'xtick',0:round(length(h)/5):length(h),'xticklabel',round(dt*(0:round(length(h)/5):length(h))))
    xlabel('time (ms)')
    title('post-spike filter')
    grid on
    
    ISI_multi = cell(1,NumSP);
    for i = 1:NumSP
        temp = zeros(NumSP,length(ISI)+NumSP-1);
        for j = 1:i
            temp(j,j:j+length(ISI)-1) = ISI;
        end
        temp = temp(:,i:end-i+1);
        ISI_multi{i} = sum(temp);
    end
    ddt = 5;
    max1 = 500;
    subplot(3,1,2);hold on
    for i = 1:NumSP
        histogram(ISI_multi{i},0:ddt:max1,'Normalization','pdf');
    end
    legend('1','2','3');
    xlim([1 length(h{i}) ]);
end


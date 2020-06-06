function dm = compileSparseDesignMatrix(dspec, trialIndices)
% Compile information from experiment according to given DesignSpec

expt = dspec.expt;
subIdxs = buildGLM.getGroupIndicesFromDesignSpec(dspec);                                        %subIdxs是1:10

totalT = sum(ceil([expt.trial(trialIndices).duration]/expt.binSize));

growingX = sparse([], [], [], 0, dspec.edim, round(totalT * dspec.edim * 0.001)); % preallocate

trialIndices = trialIndices(:)';

for kTrial = trialIndices
    nT = ceil(expt.trial(kTrial).duration / expt.binSize); % TODO move?
    
    miniX = zeros(nT, dspec.edim); % pre-allocate a dense matrix for each trial
    
    for kCov = 1:numel(dspec.covar) % for each covariate 
                                                                                                %kCov一直是1
        
        covar = dspec.covar(kCov);
        
        sidx = subIdxs{kCov};                                   %sidx一直是1:10
        
        if isfield(covar, 'cond') && ~isempty(covar.cond) && ~covar.cond(expt.trial(kTrial))
            continue;
        end
        
        stim = covar.stim(expt.trial(kTrial), nT); % either dense or sparse
        
        if isfield(covar, 'basis') && ~isempty(covar.basis)
            miniX(:, sidx) = basisFactory.convBasis(stim, covar.basis, covar.offset);

        else
            miniX(:, sidx) = stim;
        end
    end
    
    growingX = [growingX; sparse(miniX)]; %#ok<AGROW>       % growingX从0*10增加到1000*10 到2000*10。。。一直到10000*10
end

dm.X = growingX;
dm.trialIndices = trialIndices;
dm.dspec = dspec;

%% Check sanity of the design
if any(~isfinite(dm.X(:)))
    warning('Design matrix contains NaN or Inf...this is not good!');
end

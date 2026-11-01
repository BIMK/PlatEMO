function [Population,W] = MASSRMAMOEADWeights(Population,W,Z,Archive)
%MASSRMAMOEADWeights - Adapt decomposition weights using the archive.
%
%   Remove crowded subproblems, then insert sparse archive solutions
%   sequentially. New weights are normalized inverse distances from the
%   ideal point Z.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N     = numel(Population);
    M     = size(W,2);
    count = min([max(1,floor(0.05*N)),numel(Archive),N-1]);
    if count < 1
        return;
    end

    % Remove the most crowded subproblems
    for i = 1 : count
        D                    = MASSRFeatures.distance(Population.objs,Population.objs);
        D(1:size(D,1)+1:end) = inf;
        D                    = sort(D,2);
        sparsity             = prod(D(:,1:min(M,size(D,2)-1)),2);
        [~,remove]           = min(sparsity);
        Population(remove)   = [];
        W(remove,:)          = [];
    end

    % Insert sparse archive members and construct their weights
    available = true(1,numel(Archive));
    for i = 1 : count
        candidates          = find(available);
        D                   = MASSRFeatures.distance(Archive(candidates).objs,Population.objs);
        D                   = sort(D,2);
        sparsity            = prod(D(:,1:min(M,size(D,2))),2);
        [~,best]            = max(sparsity);
        selected            = candidates(best);
        Population(end+1)   = Archive(selected);
        weights             = 1./max(Archive(selected).obj-Z,1e-12);
        W(end+1,:)          = weights/sum(weights);
        available(selected) = false;
    end
end
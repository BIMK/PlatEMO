function Parents = MASSRMADEParents(Population,fitness,Problem,target,strategy)
%MASSRMADEParents - Construct the selected mating pool.
%
%   Parents contains the target individual followed by candidates from
%   random, fitness-biased, or decision-neighborhood mating, selected by
%   strategy 1, 2, or 3.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N          = numel(Population);
    candidates = setdiff(1:N,target);
    count      = min(max(6,ceil(0.15*N)),N-1);
    switch strategy
        case 1
            pool = candidates(randperm(N-1,count));
        case 2
            pool      = zeros(1,count);
            remaining = candidates;
            for j = 1 : count
                tournament                    = remaining(randperm(numel(remaining),min(3,numel(remaining))));
                [~,best]                      = min(fitness(tournament));
                pool(j)                       = tournament(best);
                remaining(remaining==pool(j)) = [];
            end
        case 3
            X                = (Population.decs-Problem.lower)./max(Problem.upper-Problem.lower,1e-12);
            distance         = sum((X-X(target,:)).^2,2);
            distance(target) = inf;
            [~,pool]         = sort(distance);
            pool             = pool(1:count)';
    end
    Parents = Population([target,pool]);
end
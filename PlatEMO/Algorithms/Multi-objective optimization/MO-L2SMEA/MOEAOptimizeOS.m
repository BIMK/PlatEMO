function candidates = MOEAOptimizeOS(subproblemList,popSize,maxIter)
% Optimize each subproblem by multi-objective algorithm

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Get dimension and boundary
    D = numel(subproblemList);
    [BU,BD] = deal(zeros(1,D));
    for k = 1 : D
        BU(k) = subproblemList{k}.ub;
        BD(k) = subproblemList{k}.lb;
    end
    
    %% Generate the random population
    population.decs = rand(popSize,D).*repmat(BU-BD,popSize,1) + repmat(BD,popSize,1);
    population.objs = SaEvaluateOS(subproblemList,population.decs);
    
    %% Optimize by multi-objective algorithm
    candidates = NSGAIIIopt(subproblemList,BU,BD,population,popSize,maxIter);
    [~,stds]   = SaEvaluateOS(subproblemList,candidates(:,1:D));
    candidates = [candidates, stds];
end
function candidates = NSGAIIIopt(subproblemList,BU,BD,Population,N,maxIter)
% Optimization via NSGA-III

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    [D,M] = deal(numel(subproblemList));
    
    %% Generate the reference points and random population
    [Z,N] = UniformPoint(N,M);
    Zmin  = min(Population.objs,[],1);

    %% Optimization
    currIter = 0;
    while currIter < maxIter
        MatingPool = TournamentSelection(2,N,sum(Population.objs,2));
        Offspring.decs = GAope(Population.decs(MatingPool,:),BU,BD);
        Offspring.objs = SaEvaluate(subproblemList,Offspring.decs);
        Zmin       = min([Zmin;Offspring.objs],[],1);
        Population = EnvironmentalSelectionIII([Population,Offspring],N,Z,Zmin); 
        currIter   = currIter + 1;
    end
    candidates = [Population.decs,Population.objs];
end
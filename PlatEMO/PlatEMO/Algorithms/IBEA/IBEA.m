function IBEA(Global)
% <algorithm> <I>
% Indicator-based evolutionary algorithm
% kappa --- 0.05 --- Fitness scaling factor

%------------------------------- Reference --------------------------------
% E. Zitzler and S. Kunzli, Indicator-based selection in multiobjective
% search, Proceedings of the International Conference on Parallel Problem
% Solving from Nature, 2004, 832-842.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    kappa = Global.ParameterSet(0.05);
    
    %% Generate random population
    Population = Global.Initialization();
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,-CalFitness(Population.objs,kappa));
        Offspring  = GA(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],Global.N,kappa);
    end
end
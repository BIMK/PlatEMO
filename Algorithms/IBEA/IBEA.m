function IBEA(Global)
% <algorithm> <H-N>
% Indicator-Based Selection in Multiobjective Search
% kappa --- 0.05 --- Fitness scaling factor

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    kappa = Global.ParameterSet(0.05);
    
    %% Generate random population
    Population = Global.Initialization();
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,-CalFitness(Population.objs,kappa));
        Offspring  = Global.Variation(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],Global.N,kappa);
    end
end
function SPEA2(Global)
% <algorithm> <O-Z>
% SPEA2: Improving the Strength Pareto Evolutionary Algorithm

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
    Fitness    = CalFitness(Population.objs);
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,Fitness);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end
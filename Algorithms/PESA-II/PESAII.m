function PESAII(Global)
% <algorithm> <O-Z>
% PESA-II: Region-based Selection in Evolutionary Multiobjective
% Optimization
% div --- 10 --- The number of divisions in each objective

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    div = Global.ParameterSet(10);

    %% Generate random population
    Population = Global.Initialization();
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,Global.N,div);
        Offspring  = Global.Variation(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],Global.N,div);
    end
end
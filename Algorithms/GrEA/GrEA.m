function GrEA(Global)
% <algorithm> <A-G>
% A Grid-Based Evolutionary Algorithm for Many-Objective Optimization
% div --- --- The number of divisions in each objective

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Div = [0 45 15 10 9 9 8 8 10 12];
    div = Global.ParameterSet(Div(min(Global.M,10)));

    %% Generate random population
    Population = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,div);
        Offspring  = Global.Variation(Population(MatingPool));    
        Population = EnvironmentalSelection([Population,Offspring],Global.N,div);
    end
end
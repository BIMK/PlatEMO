function MaOEACSS(Global)
% <algorithm> <H-N>
% Many-Objective Evolutionary Algorithms Based on Coordinated Selection
% Strategy
% t --- 0 --- Threshold value in environmental selection

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    t = Global.ParameterSet(0);

    %% Generate random population
    Population = Global.Initialization();
    Zmin       = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,Zmin);
        Offspring  = Global.Variation(Population(MatingPool));
        Zmin       = min([Zmin;Offspring.objs],[],1);
        Population = EnvironmentalSelection([Population,Offspring],Zmin,t,Global.N);
    end
end
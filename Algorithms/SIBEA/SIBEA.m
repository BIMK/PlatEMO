function SIBEA(Global)
% <algorithm> <O-Z>
% The Hypervolume Indicator Revisited: On the Design of Pareto-compliant
% Indicators Via Weighted Integration

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Liangli Zhen

    %% Generate random population
    Population = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = randi(length(Population),1,Global.N);
        Offspring  = Global.Variation(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end
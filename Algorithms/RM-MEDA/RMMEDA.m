function RMMEDA(Global)
% <algorithm> <O-Z>
% RM-MEDA: A Regularity Model-Based Multiobjective Estimation of
% Distribution Algorithm
% operator --- RMMEDA_operator

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
    
    %% Optimization
    while Global.NotTermination(Population)
        Offspring  = Global.Variation(Population,Global.N,@RMMEDA_operator);
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end
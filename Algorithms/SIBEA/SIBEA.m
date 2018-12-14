function SIBEA(Global)
% <algorithm> <S>
% Simple indicator-based evolutionary algorithm

%------------------------------- Reference --------------------------------
% E. Zitzler, D. Brockhoff, and L. Thiele, The hypervolume indicator
% revisited: On the design of Pareto-compliant indicators via weighted
% integration, Proceedings of the International Conference on Evolutionary
% Multi-Criterion Optimization, 2007, 862-876.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Liangli Zhen

    %% Generate random population
    Population = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = randi(length(Population),1,Global.N);
        Offspring  = GA(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end
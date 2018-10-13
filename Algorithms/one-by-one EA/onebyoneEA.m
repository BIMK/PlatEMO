function onebyoneEA(Global)
% <algorithm> <O-Z>
% A Many-Objective Evolutionary Algorithm Using A One-by-One Selection
% Strategy

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
    % Ideal point
    zmin = min(Population.objs,[],1);
    % Rank of each solution in one-by-one selection
    Rank = ones(1,Global.N);
    % Distribution threshold
    zeta = 1;

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,Rank);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,Rank,zeta,zmin] = EnvironmentalSelection([Population,Offspring],zeta,zmin,Global.N);
    end
end
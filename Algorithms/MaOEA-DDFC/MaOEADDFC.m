function MaOEADDFC(Global)
% <algorithm> <H-N>
% A Many-Objective Evolutionary Algorithm With Enhanced Mating and
% Environmental Selections
% K --- 5 --- The number of neighbors for estimating density
% L --- 3 --- The number of candidates for convergence-based selection

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [K,L] = Global.ParameterSet(5,3);

    %% Generate random population
    Population = Global.Initialization();
    Zmin       = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,Zmin);
        Offspring  = Global.Variation(Population(MatingPool));
        Zmin       = min([Zmin;Offspring.objs],[],1);
        Population = EnvironmentalSelection([Population,Offspring],Zmin,Global.N,K,L);
    end
end
function HypE(Global)
% <algorithm> <H-N>
% HypE: An Algorithm for Fast Hypervolume-Based Many-Objective Optimization
% nSample --- 10000 --- Number of sampled points for HV estimation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    nSample = Global.ParameterSet(10000);

    %% Generate random population
    Population = Global.Initialization();
    % Reference point for hypervolume calculation
    RefPoint = zeros(1,Global.M) + max(Population.objs)*1.2;

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,-CalHV(Population.objs,RefPoint,Global.N,nSample));
        Offspring  = Global.Variation(Population(MatingPool));    
        Population = EnvironmentalSelection([Population,Offspring],Global.N,RefPoint,nSample);
    end
end
function SMPSO(Global)
% <algorithm> <O-Z>
% SMPSO: A New PSO-based Metaheuristic for Multi-objective Optimization
% operator --- SMPSO_operator

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population       = Global.Initialization();
    Pbest            = Population;
    [Gbest,CrowdDis] = UpdateGbest(Population,Global.N);

    %% Optimization
    while Global.NotTermination(Gbest)
        Population       = Global.Variation([Population,Pbest,Gbest(TournamentSelection(2,Global.N,-CrowdDis))],Global.N,@SMPSO_operator);
        [Gbest,CrowdDis] = UpdateGbest([Gbest,Population],Global.N);
        Pbest            = UpdatePbest(Pbest,Population);
    end
end
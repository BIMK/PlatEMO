function SMPSO(Global)
% <algorithm> <S>
% Speed-constrained multi-objective particle swarm optimization

%------------------------------- Reference --------------------------------
% A. J. Nebro, J. J. Durillo, J. Garcia-Nieto, C. A. Coello Coello, F.
% Luna, and E. Alba, SMPSO: A new PSO-based metaheuristic for
% multi-objective optimization, Proceedings of the IEEE Symposium on
% Computational Intelligence in Multi-Criteria Decision-Making, 2009,
% 66-73.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population       = Global.Initialization();
    Pbest            = Population;
    [Gbest,CrowdDis] = UpdateGbest(Population,Global.N);

    %% Optimization
    while Global.NotTermination(Gbest)
        Population       = Operator(Population,Pbest,Gbest(TournamentSelection(2,Global.N,-CrowdDis)));
        [Gbest,CrowdDis] = UpdateGbest([Gbest,Population],Global.N);
        Pbest            = UpdatePbest(Pbest,Population);
    end
end
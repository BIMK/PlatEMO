classdef SMPSO < ALGORITHM
% <multi> <real/integer>
% Speed-constrained multi-objective particle swarm optimization

%------------------------------- Reference --------------------------------
% A. J. Nebro, J. J. Durillo, J. Garcia-Nieto, C. A. Coello Coello, F.
% Luna, and E. Alba, SMPSO: A new PSO-based metaheuristic for
% multi-objective optimization, Proceedings of the IEEE Symposium on
% Computational Intelligence in Multi-Criteria Decision-Making, 2009,
% 66-73.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population       = Problem.Initialization();
            Pbest            = Population;
            [Gbest,CrowdDis] = UpdateGbest(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Gbest)
                Population       = Operator(Problem,Population,Pbest,Gbest(TournamentSelection(2,Problem.N,-CrowdDis)));
                [Gbest,CrowdDis] = UpdateGbest([Gbest,Population],Problem.N);
                Pbest            = UpdatePbest(Pbest,Population);
            end
        end
    end
end
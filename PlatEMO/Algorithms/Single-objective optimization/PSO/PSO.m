classdef PSO < ALGORITHM
% <single> <real/integer> <large/none> <constrained/none>
% Particle swarm optimization
% W --- 0.4 --- Inertia weight

%------------------------------- Reference --------------------------------
% R. Eberhart and J. Kennedy, A new optimizer using particle swarm theory,
% Proceedings of the International Symposium on Micro Machine and Human
% Science, 1995, 39-43.
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
            %% Parameter setting
            W = Algorithm.ParameterSet(0.4);
            
            %% Generate random population
            Population = Problem.Initialization();
            Pbest      = Population;
            [~,best]   = min(FitnessSingle(Pbest));
            Gbest      = Pbest(best);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                Population     = OperatorPSO(Problem,Population,Pbest,Gbest,W);
                replace        = FitnessSingle(Pbest) > FitnessSingle(Population);
                Pbest(replace) = Population(replace);
                [~,best]       = min(FitnessSingle(Pbest));
                Gbest          = Pbest(best);
            end
        end
    end
end
classdef AESSPSO < ALGORITHM
% <2025> <single> <real/integer> <large/none> <constrained/none>
% Adaptive exploration state-space particle swarm optimization
% Beta  --- 2.05 --- Inertia weight
% Gamma --- 2.05 --- Inertia weight

%------------------------------- Reference --------------------------------
% M. Alimohammadi and T. M. R. Akbarzadeh. State-space adaptive exploration
% for explainable particle swarm optimization. Swarm and Evolutionary
% Computation, 2025, 94, 101868.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Mehdi Alimohammadi (alimohammadimd@chmail.ir)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [Beta,Gamma] = Algorithm.ParameterSet(2.05,2.05);

            %% Generate random population
            Population = Problem.Initialization();
            Pbest      = Population;
            [~,best]   = min(FitnessSingle(Pbest));
            Gbest      = Pbest(best);
          
            %% Optimization
            while Algorithm.NotTerminated(Population)
                Population     = OperatorAESSPSO(Problem,Population,Pbest,Gbest,Beta, Gamma);
                replace        = FitnessSingle(Pbest) > FitnessSingle(Population);
                Pbest(replace) = Population(replace);
                [~,best]       = min(FitnessSingle(Pbest));
                Gbest          = Pbest(best);
            end
        end
    end
end
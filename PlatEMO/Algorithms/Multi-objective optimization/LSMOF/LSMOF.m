classdef LSMOF < ALGORITHM
% <multi> <real/integer> <large/none>
% Large-scale multi-objective optimization framework with NSGA-II
% wD       --- 10 --- The generation of weight optimization with DE
% SubN     --- 30 --- The population size of the transferred problem
% operator ---  1 --- Original operators 1. GA 2. DE

%------------------------------- Reference --------------------------------
% C. He, L. Li, Y. Tian, X. Zhang, R. Cheng, Y. Jin, and X. Yao,
% Accelerating large-scale multi-objective optimization via problem
% reformulation, IEEE Transactions on Evolutionary Computation, 2019,
% 23(6): 949-961.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    methods
        function main(Algorithm,Problem)
            %% Parameter settings
            [wD,SubN,Operator] = Algorithm.ParameterSet(10,30,2);
            Population         = Problem.Initialization();
            G = ceil(Problem.maxFE*0.05/(SubN*2*wD));

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Problem.FE < 0.6*Problem.maxFE
                    Archive    = WeightOptimization(Problem,G,Population,wD,SubN);
                    Population = EnvironmentalSelection([Population,Archive],Problem.N);
                else
                    Population = subNSGAII(Problem,Population,Operator,Problem.N);
                end
            end
        end
    end
end
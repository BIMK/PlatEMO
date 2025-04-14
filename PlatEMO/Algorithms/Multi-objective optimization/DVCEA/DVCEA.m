classdef DVCEA < ALGORITHM
% <2025> <multi/many> <real/integer> <large/none> <constrained>
% Decision variables classification-based evolutionary algorithm

%------------------------------- Reference --------------------------------
% X. Ban, J. Liang, K. Qiao, K. Yu, Y. Wang, J. Zhu, B. Qu. A decision 
% variables classification-based evolutionary algorithm for constrained 
% multi-objective optimization problems. IEEE/CAA Journal of Automatica 
% Sinica, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population  = Problem.Initialization();
            [~, C]      = kmeans(Population.decs,5);
            [FEA,INFEA] = Variable_classification(Problem,Population,C);

            cons  = Population.cons;
            cons(cons <= 0) = 0;
            conss = sum(cons,2);
            epsilon0 = max(conss);
            if epsilon0 == 0
                epsilon0 = 1;
            end
            Fitness = CalFitness_E(Population.objs,Population.cons,epsilon0);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                cp        = (-log(epsilon0)-6)/log(1-0.5);
                epsilon   = epsilon0*(1-Problem.FE/Problem.maxFE)^cp;
                Offspring = OperatorDE_pbest_1_main(Population, Problem.N, Problem, Fitness, FEA, 0.1);
                [Population,~] = Improve_E_EnvironmentalSelection([Population,Offspring],Problem.N,epsilon);
                Offspring = DEgenerator_better(Population,Problem,INFEA,epsilon);
                [Population,Fitness] = Improve_E_EnvironmentalSelection([Population,Offspring],Problem.N,epsilon);
            end
        end
    end
end
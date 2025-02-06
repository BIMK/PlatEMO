classdef PRDH < ALGORITHM
% <2024> <multi> <binary> 
% Problem reformulation and duplication handling

%------------------------------- Reference --------------------------------
% R. Jiao, B. Xue, and M. Zhang. Solving multiobjective feature selection
% problems in classification via problem reformulation and duplication
% handling. IEEE Transactions on Evolutionary Computation, 2024, 28(4):
% 846-860.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao

    methods
        function main(Algorithm,Problem)
            %% Generate initial population
            Population = InitializePopulation(Problem);
            [~, FrontNo, CrowdDis] = EnvironmentalSelection(Population, Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);
                Offspring  = OffspringReproduction(Problem, Population(MatingPool));
                [Population, FrontNo, CrowdDis] = EnvironmentalSelection([Population, Offspring], Problem.N);
            end
        end
    end
end

function Population = InitializePopulation(Problem)
    T = min(Problem.D, Problem.N * 3);
    Pop = zeros(Problem.N, Problem.D);
    for i = 1 : Problem.N
        k = randperm(T, 1);
        j = randperm(Problem.D, k);
        Pop(i, j) = 1;
    end
    Population = Problem.Evaluation(Pop);
end
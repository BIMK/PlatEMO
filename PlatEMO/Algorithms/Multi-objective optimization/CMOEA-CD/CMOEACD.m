classdef CMOEACD < ALGORITHM
% <2025> <multi/many> <real/binary/permutation><constrained/none>
% Constraint-Pareto dominance and diversity enhancement strategy based CMOEA
% e1 --- 1 --- Type of environmental selection for forward exploration(1. SPEA2 2. NSGA-II 3. modified NSGA-III)
% e2 --- 1 --- Type of environmental selection for feasible exploitation(1. SPEA2 2. NSGA-II 3. modified NSGA-III)

%------------------------------- Reference --------------------------------
% Z. Liu, F. Han, Q. Ling, H. Han, and J. Jiang. Constraint-Pareto
% dominance and diversity enhancement strategy based evolutionary algorithm
% for solving constrained multiobjective optimization problems. IEEE
% Transactions on Evolutionary Computation, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [e1,e2] = Algorithm.ParameterSet(1,1);
            Ns = floor(Problem.N/3);

            %% Generate random population
            Population = Problem.Initialization();
            FA  = [];
            DA  = [];
            FEA = Population;
            Offspring = Population;
            zmin      = min(Population.objs,[],1) - 1e-6;

            %% Optimization
            while Algorithm.NotTerminated(FEA)
                zmin = min(zmin, min(Offspring.objs, [], 1) - 1e-6);
                FA   = ForwardExplorationArchive(FA, Offspring, zmin, Ns, e1);                   
                DA   = DiversityEnhancementArchive(DA, Offspring, zmin, Ns);
                FEA  = FeasibilityExploitationArchive(FEA, Offspring, Problem.N, e2);
                Pop1 = FA;
                Pop2 = DA;  
                Pop3 = FEA(unidrnd(length(FEA), [1, floor(Problem.N/3)]));      
                MatingPool_Pop1_1 = randperm(length(Pop1));
                MatingPool_Pop1_2 = randperm(length(Pop1));
                MatingPool_Pop2_1 = randperm(length(Pop2));
                MatingPool_Pop2_2 = randperm(length(Pop2));
                MatingPool_Pop3_1 = randperm(length(Pop3));
                MatingPool_Pop3_2 = randperm(length(Pop3));
                if rand() < 0.5
                    Offspring1 = OperatorDE(Problem, Pop1, Pop1(MatingPool_Pop1_1), Pop1(MatingPool_Pop1_2),{1,0.5,1,1});
                    Offspring2 = OperatorDE(Problem, Pop2, Pop2(MatingPool_Pop2_1), Pop2(MatingPool_Pop2_2),{1,0.5,1,1});
                    Offspring3 = OperatorDE(Problem, Pop3, Pop3(MatingPool_Pop3_1), Pop3(MatingPool_Pop3_2),{1,0.5,1,1});
                else
                    Offspring1 = OperatorGA(Problem, Pop1(MatingPool_Pop1_1), {1,20,1,1});
                    Offspring2 = OperatorGA(Problem, Pop2(MatingPool_Pop2_1), {1,20,1,1});
                    Offspring3 = OperatorGA(Problem, Pop3(MatingPool_Pop3_1), {1,20,1,1});
                end
                Offspring = [Offspring1, Offspring2, Offspring3];
            end
        end
    end
end
classdef IMTCMO < ALGORITHM
% <multi> <real/binary/permutation> <constrained>
% Improved evolutionary multitasking-based CMOEA

%------------------------------- Reference --------------------------------
% K. Qiao, J. Liang, K. Yu, C. Yue, H. Lin, D. Zhang, and B. Qu,
% Evolutionary constrained multiobjective optimization: scalable
% high-dimensional constraint benchmarks and algorithm, IEEE Transactions
% on Evolutionary Computation, 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
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
            Population1 = Problem.Initialization();
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Zmin1       = min(Population1.objs,[],1);

            Population2 = Problem.Initialization();
            Fitness2    = CalFitness(Population2.objs,Population2.cons);
            Zmin2       = min(Population2.objs,[],1);

            cons = [Population1.cons;Population2.cons];
            cons(cons<0) = 0;
            VAR0 = max(sum(cons,2));
            if VAR0 == 0
                VAR0 = 1;
            end
            X=0;
            
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                %% Udate the epsilon value
                cp  = (-log(VAR0)-6)/log(1-0.5);
                VAR = VAR0*(1-X)^cp;

                %% Offspring generation
                MatingPool = [Population1(randsample(Problem.N,Problem.N))];
                [Mate1,Mate2,Mate3]  = Neighbor_Pairing_Strategy(MatingPool,Zmin1);
                Offspring1(1:Problem.N/2) = OperatorDE_rand_1(Problem,Mate1(1:Problem.N/2), Mate2(1:Problem.N/2), Mate3(1:Problem.N/2));
                Offspring1(1+Problem.N/2:Problem.N) = OperatorDE_pbest_1_main(Population1, Problem.N/2, Problem, Fitness1, 0.1);

                MatingPool = [Population2(randsample(Problem.N,Problem.N))];
                [Mate1,Mate2,Mate3]  = Neighbor_Pairing_Strategy(MatingPool,Zmin2);
                Offspring2(1:Problem.N/2) = OperatorDE_rand_1(Problem,Mate1(1:Problem.N/2),Mate2(1:Problem.N/2),Mate3(1:Problem.N/2));
                Offspring2(1+Problem.N/2:Problem.N) = OperatorDE_pbest_1_main(Population2, Problem.N/2, Problem, Fitness2, 0.1);

                Zmin1 = min([Zmin1;Offspring1.objs],[],1);
                Zmin2 = min([Zmin2;Offspring2.objs],[],1);

                %% Environmental selection
                [Population1,Fitness1] = Main_task_EnvironmentalSelection([Population1,Offspring1,Offspring2], Problem.N, true);
                [Population2,Fitness2] = Auxiliray_task_EnvironmentalSelection([Population2,Offspring2,Offspring1], Problem.N, VAR);

                X=X+1/(Problem.maxFE/Problem.N);
            end
        end
    end
end
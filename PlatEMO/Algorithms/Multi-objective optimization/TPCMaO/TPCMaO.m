classdef TPCMaO < ALGORITHM
% <many> <real/integer/label/binary/permutation> <constrained>
% Three-population based constrained many-objective co-evolutionary algorithm

%------------------------------- Reference --------------------------------
% Y. Tian, Z. Shi, Y. Zhang, L. Zhang, H. Zhang, and X. Zhang, Solving
% optimal power flow problems via a constrained many-objective
% co-evolutionary algorithm, Frontiers in Energy Research, 2023, 11:
% 1293193.
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
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            Population3 = Problem.Initialization();
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Fitness2    = CalFitness(Population2.objs);
            Fitness3    = CalFitness(Population3.objs);

            %% Evaluate the Population
            cons = Population3.cons;
            cons(cons<0) = 0;
            epsilon0     = max(sum(cons,2));
            if epsilon0 == 0
                epsilon0 = 1;
            end
            epsilon = epsilon0;
            arch    = Archive(Population3,Problem.N);
            
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                % Evolve Population1 and Population2
                MatingPool = TournamentSelection(2,Problem.N,Fitness1);
                Offspring1 = OperatorGAhalf(Problem,Population1(MatingPool));
                MatingPool = TournamentSelection(2,Problem.N,Fitness2);
                Offspring2 = OperatorGAhalf(Problem,Population2(MatingPool));
                [tr2_1,~,feasible_rate] = EnvironmentalSelection([Population2,Offspring2],Problem.N,1);
                if feasible_rate > 0.5
                    rand_number1 = randperm(Problem.N);
                    [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,tr2_1(rand_number1(1:Problem.N/2))],Problem.N,1);
                    [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring2,Offspring1],Problem.N,2);
                else
                    [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,1);
                    [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring2,Offspring1],Problem.N,2);
                end
                % Evolve Population3
                cons = Population3.cons;
                cons(cons<0) = 0;
                if Problem.FE < 0.9*Problem.maxFE
                    if mean(sum(cons,2)<=0.02) < 0.9
                        epsilon = 0.9*epsilon;
                    else
                        epsilon = epsilon0*(1-Problem.FE/0.9/Problem.maxFE)^2;
                    end
                else
                    epsilon = 0;
                end
                MatingPool = TournamentSelection(2,Problem.N,Fitness3);
                Offspring  = OperatorGA(Problem,Population3(MatingPool));
                [Population3,Fitness3] = EnvironmentalSelection([Population3,Offspring],Problem.N,3,epsilon);
                % Output the non-dominated and feasible solutions
                arch = Archive([arch,Population3],Problem.N);
                if Problem.FE >= Problem.maxFE
                    Population1 = Archive([arch,Population1],Problem.N);
                end
            end
        end
    end
end
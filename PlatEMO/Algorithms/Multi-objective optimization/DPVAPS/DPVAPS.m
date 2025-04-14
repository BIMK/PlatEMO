classdef DPVAPS < ALGORITHM
% <2023> <multi> <real/integer> <large/none> <constrained>
% Dual-population with variable auxiliary population size
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% J. Liang, Z. Chen, Y. Wang, X. Ban, K. Qiao, and K. Yu. A dual-population
% constrained multi-objective evolutionary algorithm with variable
% auxiliary population size. Complex & Intelligent Systems, 2023, 9:
% 5907-5922.
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
            %% Parameter setting
            type = Algorithm.ParameterSet(1);

            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Fitness2    = CalFitness(Population2.objs);
            A = [];

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                if type == 1
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,size(Population2,2),Fitness2);
                    Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                    Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                elseif type == 2
                    MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                    Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                    Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                end
                cons2  = Population2.cons;
                cons2(cons2 <= 0) = 0;
                conss2 = sum(cons2,2);
                A = [A,Population2(conss2<=0)];
                if size(A,2) > Problem.N
                    [A,~] = EnvironmentalSelection(A,Problem.N,true);
                end
                [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,A,Offspring2],Problem.N,true);
                delta = (-0.9*(Problem.FE/Problem.maxFE)+1);
                [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2],round(Problem.N*delta),false);
            end
        end
    end
end
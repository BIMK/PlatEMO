classdef MTCMO < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Multitasking constrained multi-objective optimization

%------------------------------- Reference --------------------------------
% K. Qiao, K. Yu, B. Qu, J. Liang, H. Song, C. Yue, H. Lin, and K. C. Tan,
% Dynamic auxiliary task-based evolutionary multitasking for constrained
% multi-objective optimization, IEEE Transactions on Evolutionary
% Computation, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1 = Problem.Initialization();
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            
            Population2 = Problem.Initialization();
            Fitness2   = CalFitness(Population2.objs,Population2.cons);
            
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
                MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                
                MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                %% Environmental selection
                [Population1,Fitness1] = Main_task_EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,true);
                [Population2,Fitness2] = Auxiliray_task_EnvironmentalSelection([Population2,Offspring2,Offspring1],Problem.N,VAR);
                X = X + 1/(Problem.maxFE/Problem.N);
            end
        end
    end
end
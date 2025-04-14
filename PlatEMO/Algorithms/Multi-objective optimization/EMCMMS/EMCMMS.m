classdef EMCMMS < ALGORITHM
% <2024> <multi> <real/integer/label/binary/permutation> <constrained>
% Evolutionary multitasking with a cooperative multistep mutation strategy

%------------------------------- Reference --------------------------------
% K. Qiao, K. Yu, C. Yue, B. Qu, M. Liu, and J. Liang. A cooperative
% multistep mutation strategy for multiobjective optimization problems with
% deceptive constraints. IEEE Transactions on Systems, Man, and
% Cybernetics, 2024, 54(11): 6670-6682.
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
            %% Parameter settings
            run_rate    = 0.5;
            select_rate = 0.1;

            %% Initialization
            Population1 = Problem.Initialization();
            Fitness1    = CalFitness(Population1.objs,Population1.cons);

            Population2 = Problem.Initialization();
            Fitness2    = CalFitness(Population2.objs,Population2.cons);

            X    = 0;
            cons = [Population1.cons;Population2.cons];
            cons(cons<0) = 0;
            VAR0 = max(sum(cons,2));
            if VAR0 == 0
                VAR0 = 1;
            end
            cnt = 0;

            %% Optimization
            while Algorithm.NotTerminated(Population2)
                cp  = (-log(VAR0)-6)/log(1-0.5);
                VAR = VAR0*(1-X)^cp;
                cnt = cnt + 1;
                if Problem.FE/Problem.maxFE < run_rate
                    flag_index  = 2;
                    Direc_index = 2;
                    Offspring3  = CMMS( Population1, Problem, select_rate,Population2,flag_index,Direc_index);

                    [Population1,Fitness1] = Main_task_EnvironmentalSelection([Population1,Offspring3],Problem.N,true);
                    [Population2,Fitness2] = Auxiliray_task_EnvironmentalSelection( [Population2,Offspring3], Problem.N,VAR);
                end

                MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                Offspring1  = OperatorGAhalf(Problem,[Population1(MatingPool1)]);

                MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                Offspring2  = OperatorGAhalf(Problem,[Population2(MatingPool2)]);

                [Population1,Fitness1] = Main_task_EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,true);
                [Population2,Fitness2] = Auxiliray_task_EnvironmentalSelection( [Population2,Offspring2,Offspring1], Problem.N,VAR);
                X = X + 1/(Problem.maxFE/Problem.N);
            end
        end
    end
end
classdef CMOES < ALGORITHM
% <2024> <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained multi-objective optimization based on even search
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, and Y. Jin. Even search in a promising region for
% constrained multi-objective optimization. IEEE/CAA Journal of Automatica
% Sinica, 2024, 11(2): 474-486.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            type = Algorithm.ParameterSet(1);
            
            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection_NSGA2(Population1,Problem.N);
            Fitness2    = CalFitness(Population2.objs);
            Population3 = Population1;
            stage_changed = 0;
            Population = [Population1,Population2];
            CV = sum(max(0,Population.cons),2);
            max_cv = max(CV);

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                gen        = ceil(Problem.FE/Problem.N);
                CV1 = sum(max(0,Population1.cons),2);
                num_fea = sum(CV1==0);
                if num_fea <= 0 || gen <= 0.2 * ceil(Problem.maxFE/Problem.N)
                    MatingPool = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                    if type == 1
                        Offspring  = OperatorGA(Problem,Population1(MatingPool));
                    else
                        Offspring  = OperatorDE(Problem,Population1,Population1(MatingPool(1:end/2)),Population1(MatingPool(end/2+1:end)));
                    end
                    MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                    if type == 1
                        Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                    else
                        Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                    end
                    [Population1,FrontNo,CrowdDis] = EnvironmentalSelection_NSGA2([Population1,Offspring,Offspring2],Problem.N);
                    [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring,Offspring2],Problem.N,false);
                    Population3 = Population2;
                    Population = [Population1,Population2];
                    CV = sum(max(0,Population.cons),2);
                    max_cv = max(max_cv,max(CV));
                else
                    tau = gen/ceil(Problem.maxFE/Problem.N);
                    Population = [Population1,Population2,Population3];
                    CV = sum(max(0,Population.cons),2);
                    max_cv = max(CV);
                    if stage_changed == 0
                        [~,Fitness2,~,Fitness3] = EvenSearch(Population,Population1(CV1==0),Problem.N,tau,max_cv);
                        [~,Fitness1] = EnvironmentalSelection(Population1,Problem.N,true);
                        stage_changed = 1;
                    end
                    
                    if ~isempty(Population2)
                        MatingPool2 = TournamentSelection(2,2*length(Population2),Fitness2);
                        if type == 1
                            Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                        else
                            Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                        end
                    else
                        Offspring2 = [];
                    end
                    if ~isempty(Population3)
                        MatingPool3 = TournamentSelection(2,2*length(Population3),Fitness3);
                        if type == 1
                            Offspring3  = OperatorGAhalf(Problem,Population3(MatingPool3));
                        else
                            Offspring3  = OperatorDE(Problem,Population3,Population3(MatingPool3(1:end/2)),Population3(MatingPool3(end/2+1:end)));
                        end
                    else
                        Offspring3 = [];
                    end
                    MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                    if type == 1
                        Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                    else
                        Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                    end
                    
                    Offspring = [Offspring1,Offspring2,Offspring3];
                    
                    [Population1,Fitness1] = EnvironmentalSelection([Population,Offspring],Problem.N,true);
                    [Population2,Fitness2,Population3,Fitness3] = EvenSearch([Population,Offspring],Population1(CV1==0),Problem.N,tau,max_cv);
                end
            end
        end
    end
end
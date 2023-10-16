classdef CMEGL < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained evolutionary multitasking with global and local auxiliary tasks

%------------------------------- Reference --------------------------------
% K. Qiao, J. Liang, Z. Liu, K. Yu, C. Yue, and B. Qu, Evolutionary
% multitasking with global and local auxiliary tasks for constrained
% multi-objective optimization, IEEE/CAA Journal of Automatica Sinica,
% 2023, 10(10): 1951-1964.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence MaOperatorGAzine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1 = Problem.Initialization(); % Main task
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Population2 = Problem.Initialization(); % Global auxiliary task
            Fitness2   = CalFitness(Population2.objs);
            Population3 = Problem.Initialization(); % Local auxiliary task
            Fitness3   = CalFitness(Population3.objs,Population3.cons);
            
            % Calculate the constraint boundary of local auxiliary task
            cons = Population1.cons;
            cons(cons<0) = 0;
            cons =sum(cons,2);
            index =find(cons>0);
            if isempty(index)
                VAR0 = 0;
            else
                VAR0 =  mean(cons(index));
            end
            cnt  = 0;	% index of generation
            flag = 0;
            
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                cnt =cnt +1;
                if flag == 0
                    std_obj(cnt,:) = std(Population2.objs,[],1);
                    if cnt>100
                        if  sum(std(std_obj(cnt-100:cnt,:),[],1)<0.5) == Problem.M
                            flag = 1;
                        end
                    end
                end
                %% Offspring generation
                MatingPool = TournamentSelection(2,Problem.N,Fitness1);
                Offspring1 = OperatorGAhalf(Problem,[Population1(MatingPool)]);

                if flag == 0
                    MatingPool = TournamentSelection(2,Problem.N,Fitness2);
                    Offspring2 = OperatorGAhalf(Problem,[Population2(MatingPool)]);
                else
                    Offspring2 = [];
                end

                if length(Population3) <=1
                    Offspring3 = [];
                else
                    MatingPool = TournamentSelection(2,min(length(Population3),Problem.N/2),Fitness3);
                    Offspring3 = OperatorGA(Problem,[Population3(MatingPool)]);
                end

                %% Environmental selection
                [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring2,Offspring3],Problem.N,true);  
                [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1],Problem.N,true);  

                if flag == 0
                    [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2,Offspring3],Problem.N,false);
                end

                [Population3,Fitness3] = EnvironmentalSelection_LAT([Population3,Offspring1,Offspring2,Offspring3],Problem.N,VAR0);

                % Calculate the constraint boundary of local auxiliary task
                cons = Offspring1.cons;
                cons(cons<0) = 0;
                cons = sum(cons,2);
                index = find(cons>0);
                if isempty(index)
                    VAR0 = 0;
                else
                    VAR0 = mean(cons(index));
                end
            end
        end
    end
end
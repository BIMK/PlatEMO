classdef CAEAD < ALGORITHM 
% <multi> <real/integer/label/binary/permutation> <constrained>
% Dual-population evolutionary algorithm based on alternative evolution and degeneration
% type --- 1 --- Type of operator (1. DE 2. GA)

%------------------------------- Reference --------------------------------
% J. Zou, R. Sun, S. Yang, and J. Zheng, A dual-population algorithm based
% on alternative evolution and degeneration for solving constrained multi-
% objective optimization problems, Informaction Scinece, 2021, 239: 89-102.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

   methods
       function main(Algorithm,Problem)
           type = Algorithm.ParameterSet(1);
           
           %% Generate random population
           Population1 = Problem.Initialization();
           Population2 = Problem.Initialization();
           Fitness1    = CalFitness(Population1.objs,Population1.cons,0);
           Fitness2    = CalFitness(Population2.objs,Population2.cons,1e6);
           
           min_epsilon      = 1e-4;
           change_threshold = 1e-2;
           max_change       = 1;
           epsilon_k        = 1e8;
           tao              = 0.05;
           max_ep           = 0;
           gen              = 1;
           stage            = false;
           
           %% Optimization
           while Algorithm.NotTerminated(Population1)
               pop_cons2      = Population2.cons;
               cv2            = overall_cv(pop_cons2);
               population     = [Population2.decs,Population2.objs,cv2];
               Objvalues(gen) = sum(sum(Population2.objs,1));
               ep(gen)        =  epsilon_k;
               if type == 1
                   MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                   MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                   Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                   Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
               elseif type == 2
                   MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                   MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                   Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                   Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
               end
               [FrontNo2,~] = NDSort(Population2.objs,size(Population2.objs,1));
               NC2 = size(find(FrontNo2==1),2);
               if gen ~= 1
                   max_change = abs(Objvalues(gen)-Objvalues(gen-1));
               end               
               if max_change <= change_threshold &&NC2 == Problem.N && stage == false
                   epsilon_k = max(population(:,end),[],1);
                   stage = true;
               end
               Offspring3 = [];
               if stage == true
                   if type == 1
                       Offspring3 = OperatorDE(Problem,Population1,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                   elseif type == 2
                       for i=1:Problem.N/2
                           Offtemp = OperatorGAhalf(Problem,[Population1(MatingPool1(i)),Population2(MatingPool2(i))]);
                           Offspring3 = [Offspring3,Offtemp];
                       end
                   end
               end
               if stage == true
                   [stage,epsilon_k] =  update_epsilon(stage,tao,epsilon_k,max_ep,min_epsilon);
               end
               if epsilon_k < 9e5
                   max_ep = max(max_ep,epsilon_k);
               end
               [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2,Offspring3],Problem.N,true,0);
               [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring2],Problem.N,false,epsilon_k);
               gen = gen+1;
           end
       end
   end
end

function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end

function [stage,result] = update_epsilon(stage,tao,epsilon_k,epsilon_0,min_epsilon)
    if epsilon_k > min_epsilon
        result = (1 - tao) * epsilon_k;
        stage = true;
    else
        result = epsilon_0;
        stage = false;
    end
end
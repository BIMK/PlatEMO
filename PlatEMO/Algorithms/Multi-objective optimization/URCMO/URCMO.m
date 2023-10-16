classdef URCMO < ALGORITHM
% <multi> <real/integer> <constrained>
% Utilizing the relationship between constrained and unconstrained Pareto fronts for constrained multi-objective optimization

%------------------------------- Reference --------------------------------
% J. Liang, K. Qiao, K. Yu, B. Qu, C. Yue, W. Guo, and L. Wang, Utilizing
% the relationship between unconstrained and constrained Pareto fronts for
% constrained multi-objective optimization, IEEE Transactions on
% Cybernetics, 2022.
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
            warning off
            Population{1} = Problem.Initialization();
            Population{2} = Problem.Initialization();
            
            Fitness{1} = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{2} = CalFitness(Population{2}.objs);
            
            %% parameter settings
            cnt       = 0;
            p         = 0.1;
            first_FES = 10000;
            beita     = 0.9;
            while Algorithm.NotTerminated(Population{1})
                cnt = cnt + 1;
                %% learning phase
                if Problem.FE < first_FES
                    for i = 1 : 2
                        valOffspring{i}(1:Problem.N/2)             = GA_TournamentSelection(Problem, Population{i}, Fitness{i}, Problem.N/2);
                        valOffspring{i}(1+Problem.N/2 : Problem.N) = DE_current_to_rand_1(Problem, Population{i}, Problem.N/2);
                    end
                    
                    % Environmental selection
                    [Population{1},succ1,succ2,Fitness{1}] = First_Stage_EnvironmentalSelection([Population{1},valOffspring{1},valOffspring{2}],Problem.N,1);
                    [Population{2},~,~,Fitness{2}]         = First_Stage_EnvironmentalSelection([Population{2},valOffspring{2},valOffspring{1}],Problem.N,2);
                    
                    succ_jilu(cnt,1) = sum(succ1(1:Problem.N/2));
                    succ_jilu(cnt,2) = sum(succ1(1 + Problem.N/2 : Problem.N));
                    
                    succ_jilu(cnt,3) = sum(succ2(1:Problem.N/2));
                    succ_jilu(cnt,4) = sum(succ2(1 + Problem.N/2 : Problem.N)); % only succ2 is used in the paper
                    
                    if Problem.FE >= first_FES
                        for num = 1 : 4
                            if std(succ_jilu(:,num)) ~= 0
                                a(num) = mean(succ_jilu(:,num)) / std(succ_jilu(:,num));
                            else
                                a(num) = 0;
                            end
                        end
                        
                        [flag,ll] = Classification(Population{1},Population{2},beita);  % flag correspondings to the three kinds of evolving strategies
                        
                        if flag == 1
                            if a(3) < a(4)
                                flag = 3;
                            end
                        end
                        
                    end
                    
                    %% Evolving phase
                else
                    
                    valOffspring{1}(1:Problem.N/2)             = GA_TournamentSelection(Problem, Population{1}, Fitness{1}, Problem.N/2);
                    valOffspring{1}(1+Problem.N/2 : Problem.N) = DE_current_to_rand_1(Problem, Population{1}, Problem.N/2);
                    
                    if flag == 1
                        valOffspring{2}(1:Problem.N/2)              = GA_TournamentSelection(Problem, Population{2}, Fitness{2}, Problem.N/2);
                        valOffspring{2}( 1+Problem.N/2 : Problem.N) = DE_transfer(Problem, Population{2}, Population{1}, Problem.N/2);
                    elseif flag == 2
                        num      = ceil(Problem.N/3);
                        rand_num = rand;
                        if rand_num <= 1/3
                            valOffspring{2}(1:num)               = GA_TournamentSelection(Problem, Population{2}, Fitness{2}, num);
                            valOffspring{2}(1+num : 2*num)       = DE_transfer(Problem, Population{2}, Population{1}, num);
                            valOffspring{2}(1+2*num : Problem.N) = DE_current_to_other_pbest_1(Problem,Population{2},Problem.N-2*num,Fitness{1},Population{1},p);
                        elseif rand_num <=2/3
                            valOffspring{2}(1:num)               = GA_TournamentSelection(Problem, Population{2}, Fitness{2}, num);
                            valOffspring{2}(1+num : 2*num)       = DE_current_to_other_pbest_1(Problem,Population{2},num,Fitness{1},Population{1},p);
                            valOffspring{2}(1+2*num : Problem.N) = DE_transfer(Problem, Population{2}, Population{1}, Problem.N-2*num);
                        elseif rand_num <= 1
                            valOffspring{2}(1:num)               = DE_current_to_other_pbest_1(Problem,Population{2},num,Fitness{1},Population{1},p);
                            valOffspring{2}(1+num : 2*num)       = DE_transfer(Problem, Population{2}, Population{1}, num);
                            valOffspring{2}(1+2*num : Problem.N) = GA_TournamentSelection(Problem, Population{2}, Fitness{2}, Problem.N-2*num);
                        end
                    elseif flag == 3
                        valOffspring{2} = DE_current_to_other_pbest_1(Problem,Population{2},Problem.N,Fitness{1},Population{1},p);
                    end
                    
                    % Environmental selection
                    [Population{1},Fitness{1}] = Second_Stage_EnvironmentalSelection([Population{1},valOffspring{1},valOffspring{2}],Problem.N,1);
                    [Population{2},Fitness{2}] = Second_Stage_EnvironmentalSelection([Population{2},valOffspring{2},valOffspring{1}],Problem.N,2);
                end
            end
        end
    end
end
classdef EMCMO < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Evolutionary multitasking-based constrained multiobjective optimization

%------------------------------- Reference --------------------------------
% K. Qiao, K. Yu, B. Qu, J. Liang, H. Song, and C. Yue. An evolutionary
% multitasking optimization framework for constrained multi-objective
% optimization problems, IEEE Transactions on Evolutionary Computation,
% 2022, 26(2): 263-277.
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
            Population{1} = Problem.Initialization();
            Population{2} = Problem.Initialization();
            Fitness{1}    = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{2}    = CalFitness(Population{2}.objs);
            transfer_state=0;
            
            cnt=0;
            %% Optimization
            while Algorithm.NotTerminated(Population{1})
                cnt =cnt+1;
                if   transfer_state == 0
                    for i = 1: 2
                        valOffspring{i} = OperatorGAhalf(Problem,Population{i}(randi(Problem.N,1,Problem.N)));
                    end
                    
                    for i = 1:2
                        if i == 1
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{1},valOffspring{1},valOffspring{2}],Problem.N,i);
                        else
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{2},valOffspring{2},valOffspring{1}],Problem.N,i);
                        end
                    end
                    
                    if Problem.FE/Problem.maxFE >=0.2
                        transfer_state = 1;
                    end
                    
                else
                    
                    for i = 1: 2
                        MatingPool = TournamentSelection(2,Problem.N,Fitness{i});
                        valOffspring{i} = OperatorGAhalf(Problem,Population{i}(MatingPool));
                    end
                    [~,~,Next] = EnvironmentalSelection( [Population{2},valOffspring{2}],Problem.N,1);
                    succ_rate(1,cnt) =  (sum(Next(1:Problem.N))/100) - (sum(Next(Problem.N+1:end))/50);
                    
                    [~,~,Next] = EnvironmentalSelection( [Population{1},valOffspring{1}],Problem.N,2);
                    succ_rate(2,cnt) =  (sum(Next(1:Problem.N))/100) - (sum(Next(Problem.N+1:end))/50);
                    
                    for i = 1:2
                        if   succ_rate(i,cnt) >0
                            rand_number = randperm(Problem.N);
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{i},valOffspring{i},Population{2/i}(rand_number(1:Problem.N/2))],Problem.N,i);
                        else
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{i},valOffspring{i},valOffspring{2/i}],Problem.N,i);
                        end
                    end
                end
            end
        end
    end
end
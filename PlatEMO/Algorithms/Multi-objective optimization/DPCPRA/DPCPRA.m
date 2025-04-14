classdef DPCPRA < ALGORITHM
% <2024> <multi> <real/integer/label/binary/permutation> <constrained>
% Dual-population with dynamic constraint processing and resource allocating

%------------------------------- Reference --------------------------------
% K. Qiao, Z. Chen, B. Qu, K. Yu, C. Yue, K. Chen, and J. Liang. A dual-
% population evolutionary algorithm based on dynamic constraint processing
% and resources allocation for constrained multi-objective optimization
% problems. Expert Systems With Applications, 2024, 238: 121707.
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
            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            Fitness1    = CalFitness_pop1(Population1.objs,Population1.cons);
            current_cons       = 0;
            gen                = 0;
            last_gen           = 100;
            change_threshold   = 1e-2;
            change_rate        = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            priority           = [];
            flag               = 0;
            constraint_handing = 0;
            archive            = Population2;
            Fitness2           = CalFitness_pop2(Population2.objs,Population2.cons,priority,current_cons,constraint_handing);
            success_rate1      = 0.5;

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                if flag == 0
                    change_rate = Normalization(Population2,change_rate,ceil(Problem.FE/Problem.N));
                    if Convertion(change_rate,ceil(Problem.FE/Problem.N),gen,last_gen,change_threshold)
                        flag = 1;
                        [priority,evaluatedasible_rate] = Constraint_priority(Population2);
                        Population2 = Problem.Initialization();
                    end
                else
                    % Judge whether to enter next stage
                    if current_cons == 0
                        CV = Population2.cons;
                        CV = CV(:,priority(1));
                        if length(find(CV>0))/Problem.N > 0
                            current_cons = current_cons + 1;
                            gen = ceil(Problem.FE/Problem.N) + 1;
                        end
                    elseif current_cons <= size(Population2.cons,2)
                        if constraint_handing == 0
                            change_rate = Normalization(Population2,change_rate,ceil(Problem.FE/Problem.N));
                            if Convertion(change_rate,ceil(Problem.FE/Problem.N),gen,last_gen,change_threshold)
                                if current_cons<size(Population2.cons,2) && evaluatedasible_rate(priority(current_cons+1))~=1
                                    current_cons = current_cons+1;
                                elseif current_cons<size(Population2.cons,2) && evaluatedasible_rate(priority(current_cons+1))==1
                                    current_cons = size(Population2.cons,2);
                                elseif current_cons == size(Population2.cons,2)
                                    constraint_handing = 1;
                                end
                                if size(archive,2) == Problem.N
                                    Population2 = archive;
                                    Fitness2    = CalFitness_pop2(Population2.objs,Population2.cons,priority,current_cons,constraint_handing);
                                else
                                    archive     = Archive([archive,Population2],Problem.N,priority,current_cons,size(archive,2));
                                    Population2 = archive;
                                    Fitness2    = CalFitness_pop2(Population2.objs,Population2.cons,priority,current_cons,constraint_handing);
                                end
                                gen = ceil(Problem.FE/Problem.N) + 1;
                            end
                        end
                    end
                end
                % Optimization
                if flag == 0
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                    Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                    Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                else
                    if mod(ceil(Problem.N*success_rate1),2)~=0
                        MatingPool1 = TournamentSelection(2,ceil(Problem.N*success_rate1)+1,Fitness1);
                    else
                        MatingPool1 = TournamentSelection(2,ceil(Problem.N*success_rate1),Fitness1);
                    end

                    Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                    MatingPool2 = TournamentSelection(2,Problem.N-2*length(Offspring1),Fitness2);
                    Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                end
                %  Update external archive
                if flag==1 && constraint_handing~=1
                    archive = Archive([Offspring2,archive],Problem.N,priority,current_cons);
                end
                [Population1,Fitness1,success_rate1] = EnvironmentalSelection_pop1([Population1,Offspring1,Offspring2],Problem.N,true,length(Offspring1));
                [Population2,Fitness2,success_rate2] = EnvironmentalSelection_pop2([Population2,Offspring2,Offspring1],Problem.N,priority,current_cons,constraint_handing,length(Offspring2));
                success_rate1 = success_rate1/(success_rate1+success_rate2);
            end
        end
    end
end
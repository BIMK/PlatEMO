classdef MSCMO < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Multi-stage constrained multi-objective evolutionary algorithm
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% H. Ma, H. Wei, Y. Tian, R. Cheng, and X. Zhang, A multi-stage
% evolutionary algorithm for multi-objective optimization with complex
% constraints, Information Sciences, 2021, 560: 68-91.
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
            %% Parameter setting
            type = Algorithm.ParameterSet(1);

           %% Initialization
            Population         = Problem.Initialization();
            gen                = 0;
            last_gen           = 100;
            change_threshold   = 1e-2;
            change_rate        = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            priority           = [];
            current_cons       = 0;
            flag               = 0;
            constraint_handing = 0;
            Fitness            = CalFitness(Population.objs,Population.cons,priority,current_cons,constraint_handing);
            archive            = Population;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Determine constraint-handling priority
                if flag == 0
                    change_rate = Normalization(Population,change_rate,ceil(Problem.FE/Problem.N));
                     if Convertion(change_rate,ceil(Problem.FE/Problem.N),gen,last_gen,change_threshold)
                        flag = 1;
                        [priority,Feasible_rate] = Constraint_priority(Population);
                        Population = Problem.Initialization();
                     end
                else
                    % Judge whether to enter next stage
                    if current_cons == 0
                       CV = Population.cons;
                       CV = CV(:,priority(1));
                       if length(find(CV>0))/Problem.N > 0
                          current_cons = current_cons + 1;
                          gen = ceil(Problem.FE/Problem.N) + 1;
                       end
                    elseif current_cons <= size(Population.cons,2)
                         if constraint_handing == 0
                             change_rate = Normalization(Population,change_rate,ceil(Problem.FE/Problem.N));
                             if Convertion(change_rate,ceil(Problem.FE/Problem.N),gen,last_gen,change_threshold)
                                 if current_cons<size(Population.cons,2) && Feasible_rate(priority(current_cons+1))~=1
                                    current_cons = current_cons+1;
                                 elseif current_cons<size(Population.cons,2) && Feasible_rate(priority(current_cons+1))==1
                                    current_cons = size(Population.cons,2);
                                 elseif current_cons == size(Population.cons,2)
                                    constraint_handing = 1;
                                 end
                                 if size(archive,2) == Problem.N
                                    Population = archive;
                                    Fitness    = CalFitness(Population.objs,Population.cons,priority,current_cons,constraint_handing);
                                 else
                                    archive    = Archive([archive,Population],Problem.N,priority,current_cons,size(archive,2));
                                    Population = archive;
                                    Fitness    = CalFitness(Population.objs,Population.cons,priority,current_cons,constraint_handing);
                                 end
                                 gen = ceil(Problem.FE/Problem.N) + 1;
                             end
                         end  
                    end     
                end
                % Reproduction
                if type == 1
                    MatingPool = TournamentSelection(2,Problem.N,Fitness);
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                elseif type == 2
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness);
                    MatingPool2 = TournamentSelection(2,Problem.N,Fitness);
                    Offspring   = OperatorDE(Problem,Population,Population(MatingPool1),Population(MatingPool2));
                end
                %  Update external archive
                if flag==1 && constraint_handing~=1
                    archive = Archive([Offspring,archive],Problem.N,priority,current_cons);
                end
                % Environmental selection 
                [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Problem.N,priority,current_cons,constraint_handing);
            end
        end
    end
end
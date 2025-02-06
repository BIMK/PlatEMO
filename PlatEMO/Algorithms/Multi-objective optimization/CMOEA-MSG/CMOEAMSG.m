classdef CMOEAMSG < ALGORITHM
% <2024> <multi> <real/integer> <constrained>
% Multi-stage constrained multi-objective evolutionary algorithm
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% Y. Tian, J. Chen, and X. Zhang. An optimizer combining evolutionary
% computation and gradient descent for constrained multi-objective
% optimization. Journal of Computer Applications (Chinese), 2024, 44(05):
% 1386-1392.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
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
            priority     = [];
            current_cons = 0;
            Population1  = Problem.Initialization(100);
            Population2  = Problem.Initialization(100);
            Population   = [Population1,Population2];
            Ar = Population1;       
            i  = 0;
            ii = 0;
            ab = 10;
            gen                = 0;
            change_threshold   = 0.1;
            change_rate        = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            constraint_handing = 0;
            [priority,~]       = Constraint_priority(Population,priority,current_cons);
            last_gen           = max(30,floor(200/(length(priority)+1)));
            
            %% Optimization
            Population1 = EnvironmentalSelection1(Population1,100,priority,current_cons,0);
            for h = 1 : length(Population1)
                while true
                    X    = Population1(h);
                    Grad = Gradient(Problem,X,current_cons,priority,0);
                    Grad = sign(Grad);
                    while true
                        step = abs(Population1(h).dec-Population1(randi(100,1)).dec);
                        i = i + 1;
                        H = X;
                        X = X.dec - step.*Grad*0.1  ;
                        X = Problem.Evaluation(X);
                        [X,Bool] = Compare([X,H],priority,current_cons,constraint_handing);
                        if i ==ab || Bool == false
                            break;
                        end
                    end
                    if i == 10 || Bool == false
                        i = 0;
                        break;
                    end
                end
                Population1 = [Population1,X];
            end
            [Population1,~] = EnvironmentalSelection(Population1,100,priority,current_cons,0);
            while Algorithm.NotTerminated(Ar)
                ii = ii + 1;
                if current_cons <= size(Population1.cons,2)
                    if constraint_handing == 0
                        change_rate(ii,:) = mean(Ar.objs);
                        if (ii-gen > last_gen) && (min(abs(change_rate(ii,:)-mean(change_rate(ii-30:ii-1,:))))<=change_threshold)
                            if current_cons < size(Population1.cons,2)
                                [priority,~] = Constraint_priority(Population1,priority,current_cons);
                                current_cons = current_cons+1;
                            elseif current_cons == size(Population1.cons,2)
                                constraint_handing = 1;
                            end
                            gen = ii + 1;
                        end
                    end
                end
                if constraint_handing == 0
                    if  mod(ii,1) == 0
                        Fitness1 = CalFitness(Population1.objs,Population1.cons,priority,current_cons,constraint_handing);
                        [Xx,oc]  = Selection(Population1,current_cons,Fitness1,priority);
                        if current_cons == 0
                            m = 3;
                        else
                            m = 1;
                        end
                        for h = 1 : m
                            X    = Xx(h);
                            Grad = Gradient(Problem,X,current_cons,priority,oc);
                            Grad = sign(Grad);
                            while true
                                i    = i + 1;
                                step = abs(Population1(h).dec-Population1(randi(100,1)).dec);
                                H    = X;
                                X    = X.dec - step.*Grad*0.1  ;
                                X    = Problem.Evaluation(X);
                                [X,Bool] = Compare([X,H],priority,current_cons,constraint_handing);
                                if i > ab || Bool == false
                                    i = 0;
                                    break;
                                end
                            end
                            Population1 = [Population1,X];
                        end
                    end
                end
                Fitness = CalFitness(Population1.objs,Population1.cons,priority,length(priority),0);
                if type == 1
                    MatingPool1 = TournamentSelection(2,100,Fitness);
                    Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                else
                    a = length(Population1);
                    MatingPool11 = TournamentSelection(2,a,Fitness);
                    MatingPool12 = TournamentSelection(2,a,Fitness);
                    Offspring1   = OperatorDE(Problem,Population1,Population1(MatingPool11),Population1(MatingPool12));
                end
                Fitness5 = CalFitness(Ar.objs,Ar.cons,priority,length(priority),1);
                if type == 1
                    MatingPool4 = TournamentSelection(2,Problem.N,Fitness5);
                    Offspring4  = OperatorGA(Problem,Ar(MatingPool4));
                    if constraint_handing == 0
                        [Population1,~] = EnvironmentalSelection([Population1,Offspring1],100,priority,current_cons,0);
                    else
                        [Population1,~] = EnvironmentalSelection1([Population1,Offspring1],100,priority,current_cons,0);
                    end
                    Ar = EnvironmentalSelection1([Ar,Offspring4,Offspring1],Problem.N,priority,length(priority),1);
                else
                    MatingPool41 = TournamentSelection(2,100,Fitness5);
                    MatingPool42 = TournamentSelection(2,100,Fitness5);
                    Offspring4   = OperatorDE(Problem,Ar,Ar(MatingPool41),Ar(MatingPool42));
                    if current_cons == 0
                        [Population1,~] = EnvironmentalSelection([Population1,Offspring1],100,priority,current_cons,0);
                    else
                        [Population1,~] = EnvironmentalSelection1([Population1,Offspring1],100,priority,current_cons,0);
                    end
                    Ar = EnvironmentalSelection1([Ar,Offspring4,Offspring1],Problem.N,priority,length(priority),1);
                end
            end
        end
    end
end
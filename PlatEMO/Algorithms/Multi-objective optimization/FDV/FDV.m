classdef FDV < ALGORITHM
% <multi/many> <real/integer> <large/none> 
% Fuzzy decision variable framework with various internal optimizers
% Rate      --- 0.8 --- Fuzzy evolution rate. Default = 0.8
% Acc       --- 0.4 --- Step acceleration. Default = 0.4
% optimizer ---   5 --- Internal optimisation algorithm. 1 = NSGA-II, 2 = NSGA-III, 3 = MOEA/D, 4 = CMOPSO, 5 = LMOCSO
% type      ---   1 --- The type of aggregation function for MOEA/D

%------------------------------- Reference --------------------------------
% X. Yang, J. Zou, S. Yang, J. Zheng, and Y. Liu, A fuzzy decision
% variables framework for large-scale multiobjective optimization, IEEE
% Transactions on Evolutionary Computation, 2021.
%--------------------------------------------------------------------------

%  Copyright (C) 2021 Xu Yang
%  Xu Yang <xuyang.busyxu@qq.com> or <xuyang369369@gmail.com>

    methods
        function main(Algorithm,Problem)
            %% Set the default parameters
            [Rate,Acc,optimizer,type] = Algorithm.ParameterSet(0.8,0.4,5,1);

            %% NSGAII
            if optimizer==1
                % Generate random population
                Population = Problem.Initialization();
                [~,FrontNo,CrowdDis] = EnvironmentalSelection_NSGAII(Population,Problem.N);

                % Optimization
                while Algorithm.NotTerminated(Population)
                    MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                    OffDec     = OperatorGA(Problem,Population(MatingPool).decs);
                    %% FDV
                    if Problem.FE/Problem.maxFE <= Rate
                        Offspring = FDVOperator(Problem,Rate,Acc,OffDec);
                    else
                        Offspring = Problem.Evaluation(OffDec);
                    end
                    %% 
                    [Population,FrontNo,CrowdDis] = EnvironmentalSelection_NSGAII([Population,Offspring],Problem.N);
                end
            end

            %% NSGAIII
            if optimizer==2
                % Generate the reference points and random population
                [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
                Population    = Problem.Initialization();
                Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);

                % Optimization
                while Algorithm.NotTerminated(Population)
                    MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                    OffDec     = OperatorGA(Problem,Population(MatingPool).decs);
                    %% FDV
                    if Problem.FE/Problem.maxFE <= Rate
                        Offspring = FDVOperator(Problem,Rate,Acc,OffDec);
                    else
                        Offspring = Problem.Evaluation(OffDec);
                    end
                    Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                    Population = EnvironmentalSelection_NSGAIII([Population,Offspring],Problem.N,Z,Zmin);
                end
            end

            %% MOEA/D
            if optimizer==3
                % Generate the weight vectors
                [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
                T = ceil(Problem.N/10);

                % Detect the neighbours of each solution
                B = pdist2(W,W);
                [~,B] = sort(B,2);
                B = B(:,1:T);

                % Generate random population
                Population = Problem.Initialization();
                Z = min(Population.objs,[],1);

                % Optimization
                while Algorithm.NotTerminated(Population)
                    % For each solution
                    for i = 1 : Problem.N      
                        % Choose the parents
                        P = B(i,randperm(size(B,2)));

                        % Generate an offspring
                        OffDec = OperatorGAhalf(Problem,Population(P(1:2)).decs);
                        %% FDV
                        if Problem.FE/Problem.maxFE <= Rate
                            Offspring = FDVOperator(Problem,Rate,Acc,OffDec);
                        else
                            Offspring = Problem.Evaluation(OffDec);
                        end

                        % Update the ideal point
                        Z = min(Z,Offspring.obj);

                        % Update the neighbours
                        switch type
                            case 1
                                % PBI approach
                                normW   = sqrt(sum(W(P,:).^2,2));
                                normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
                                normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                                CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                                CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
                                g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                                g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                            case 2
                                % Tchebycheff approach
                                g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                                g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
                            case 3
                                % Tchebycheff approach with normalization
                                Zmax  = max(Population.objs,[],1);
                                g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                                g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                            case 4
                                % Modified Tchebycheff approach
                                g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                                g_new = max(repmat(abs(Offspring.obj-Z),T,1)./W(P,:),[],2);
                        end
                        Population(P(g_old>=g_new)) = Offspring;
                    end
                end
            end

            %% CMOPSO
            if optimizer == 4
                % Generate random population
                Population = Problem.Initialization();

                % Optimization
                while Algorithm.NotTerminated(Population)
                    [OffDec,OffVel] = Operator_CMOPSO(Problem,Population);
                    %% FDV
                    if Problem.FE/Problem.maxFE <= Rate
                        Offspring = FDVOperator(Rate,Acc,OffDec,OffVel);
                    else
                        Offspring = Problem.Evaluation(OffDec,OffVel);
                    end
                    Population = EnvironmentalSelection_CMOPSO([Population,Offspring],Problem.N);
                end
            end

            %% LMOCSO
            if optimizer == 5
                 % Generate random population
                [V,Problem.N] = UniformPoint(Problem.N,Problem.M);
                Population    = Problem.Initialization();
                Population    = EnvironmentalSelection_LMOCSO(Population,V,(Problem.FE/Problem.maxFE)^2);

                % Optimization
                while Algorithm.NotTerminated(Population)
                    % Calculate the fitness by shift-based density   SDE (the shift-based density estimation strategy)
                    PopObj = Population.objs;
                    N      = size(PopObj,1);
                    fmax   = max(PopObj,[],1);
                    fmin   = min(PopObj,[],1);
                    PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
                    Dis    = inf(N);
                    for i = 1 : N
                        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
                        for j = [1:i-1,i+1:N]
                            Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:)); 
                        end
                    end
                    Fitness = min(Dis,[],2); 

                    if length(Population) >= 2
                        Rank = randperm(length(Population),floor(length(Population)/2)*2);
                    else
                        Rank = [1,1];
                    end
                    Loser  = Rank(1:end/2);
                    Winner = Rank(end/2+1:end);
                    Change = Fitness(Loser) >= Fitness(Winner);
                    Temp   = Winner(Change);
                    Winner(Change) = Loser(Change);
                    Loser(Change)  = Temp;

                    [OffDec,OffVel] = Operator_LMOCSO(Problem,Population(Loser),Population(Winner),Rate);
                    %% FDV
                    iter = Problem.FE/Problem.maxFE;
                    if iter <= Rate
                        Offspring = FDVOperator(Problem,Rate,Acc,OffDec,OffVel);
                    else
                        Offspring = Problem.Evaluation(OffDec,OffVel);
                    end
                    Population = EnvironmentalSelection_LMOCSO([Population,Offspring],V,(Problem.FE/Problem.maxFE)^2);
                end
            end
        end
    end
end
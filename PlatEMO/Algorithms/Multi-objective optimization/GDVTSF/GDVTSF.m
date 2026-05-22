classdef GDVTSF < ALGORITHM
% <2025> <multi> <real> <large/none>
% Generational difference vector based tri-entropy structure framework
% delta --- 0.3 --- The threshold ratio for WOF framework
% K     ---   5 --- The number of clusters
% B     ---  20 --- The number of sampling points
% G     --- 0.5 --- Global control factor
% T_clu ---  20 --- Max iterations for clustering

%------------------------------- Reference --------------------------------
% Y. Xu, Y. Zhang, and W. Hu. A generational difference vector based
% tri-entropy structure optimizer for large-scale multiobjective
% optimization. Swarm and Evolutionary Computation, 2025, 98: 102079.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Hannyx Xu
% If you have any questions, please open an issue in our GitHub repository:
% https://github.com/yuzhang576/GDVTSF

    methods
        function main(Algorithm,Problem)
            %% Core Parameters Setting (Analyzed in the paper)
            [delta, K, B, G, T_clu] = Algorithm.ParameterSet(0.3, 5, 20, 0.5, 20);
            
            %% Default Settings for WOF and FDV
            [Rate, Acc]          = deal(0.6, 0.4);
            [gamma, groups, psi] = deal(4, 2, 3);
            [t1, t2]             = deal(1000, 500);
            q                    = Problem.M + 1;
            t_pop                = 10;
            sel_m                = 3;

            %% Generate a set of uniformly distributed points
            [V, Problem.N] = UniformPoint(Problem.N, Problem.M);
            N              = Problem.N;
            D              = Problem.D;
            
            %% Initialization
            % Generate random population
            population = Problem.Initialization();
            LpopObj    = population.objs;
            LpopDec    = population.decs;
            swarmN     = floor(N/3);
            % Global best solution
            gBest = GDVTSF_Utils.update_gbest(population, swarmN, V, 0);
            % Clustering decision vectors
            [Lidx, Lcenter] = kmeans(LpopDec, K, "MaxIter", T_clu);
            Lrd = zeros(K,D);
            for i = 1 : K
                ch       = Lidx == i;
                Lrd(i,:) = mean(abs(Lcenter(i,:)-LpopDec(ch,:)),1);
            end
            % Generational difference vector
            GDV = zeros(N,D);
            % Depart to 3 swarms randomly
            if length(population) >= 3
                Rank = randperm(length(population),swarmN*3);
            else
                Rank = [1,1,1];
            end
            p1 = Rank(1:swarmN);
            p2 = Rank(swarmN+1:2*swarmN);
            p3 = Rank(2*swarmN+1:end);
            % Initialize WOF Environment
            WOF_WeightIndividual.Current(Problem);
            
            %% Optimization loop
            while Algorithm.NotTerminated(gBest)
                if Problem.FE < delta*Problem.maxFE
                    population = WOF_Utils.fill_population(population, Problem);
                    % Normal optimisation step for t1 evaluations
                    population = GDVTSF_Utils.WOF_optimizer_GDV_TSO(Problem, population, t1, K, B, G, T_clu); 
                    % Selection of xPrime solutions 
                    xPrimeList = WOF_Utils.select_xprimes(population, q, sel_m); 
                    WList      = [];
                    % Do for each xPrime
                    for c = 1 : size(xPrimeList,2)
                        xPrime = xPrimeList(c);
                        % Create variable groups 
                        [Group, gamma_n] = WOF_Utils.create_groups(Problem, gamma, xPrime, groups);
                        % A dummy object is needed to simulate the global class
                        GlobalDummy = WOF_Utils.create_global_dummy(gamma_n, xPrime, Group, Problem, t_pop, psi, 2);
                        % Create initial population for the transformed problem
                        WeightPopulation = WOF_Utils.create_initial_weight_population(GlobalDummy.N, gamma_n, GlobalDummy);
                        % Optimise the transformed problem 
                        WeightPopulation = WOF_Utils.optimise_by_MOEAD(GlobalDummy, WeightPopulation, GlobalDummy.uniW, t2-t_pop, true);
                        % Extract the population 
                        W     = WOF_Utils.extract_population(WeightPopulation, Problem, population, Group, psi, xPrime, q, sel_m);
                        WList = [WList,W];
                    end
                    % Join populations. Duplicate solutions need to be removed. 
                    population = WOF_Utils.eliminate_duplicates([population,WList]);
                    population = WOF_Utils.fill_population(population, Problem);
                    % Environmental Selection
                    [population,~,~] = WOF_Utils.environmental_selection(population,Problem.N);
                else
                    % Generate new population
                    popDec = GDVTSF_Utils.operator_GDV_TSO(Problem, population, gBest, GDV, p1, p2, p3, B, G);
                    popObj = population.objs;
                    % Calculate GDV
                    [Lcenter, Lidx, Lrd, GDV] = GDVTSF_Utils.cal_GDV(N, D, K, T_clu, popDec, popObj, Lrd, LpopObj, Lcenter, Lidx);
                    LpopObj                   = popObj;
                    % FDV (Fuzzy Difference Vector operation)
                    iter = Problem.FE/Problem.maxFE;
                    if iter <= Rate
                        population = GDVTSF_Utils.operator_FDV(Problem,Rate,Acc,popDec);
                    else
                        population = Problem.Evaluation(popDec);
                    end
                end
                % Update gBest
                gBest = GDVTSF_Utils.update_gbest([gBest, population], swarmN, V, (Problem.FE/Problem.maxFE)^2);
            end
        end
    end
end
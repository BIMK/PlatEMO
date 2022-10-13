classdef WOF < ALGORITHM
% <multi> <real/integer> <large/none>
% Weighted optimization framework
% gamma         --- 4    --- Number of groups. Default = 4 
% groups        --- 2    --- Grouping method, 1 = linear, 2 = ordered, 3 = random. Default = ordered 
% psi           --- 3    --- Transformation function, 1 = Multiplication, 2 = P-Value, 3 = Interval. Default = Interval
% t1            --- 1000 --- Number of evaluations for original problem. Default = 1000
% t2            --- 500  --- Number of evaluations for transformed problem. Default = 500
% q             ---      --- The number of chosen solutions to do weight optimisation. If no value is specified, the default value is M+1
% delta         --- 0.5  --- The fraction of function evaluations to use for the alternating weight-optimisation phase. Default = 0.5
% optimiser     --- 1    --- Internal optimisation algorithm. 1 = SMPSO, 2 = MOEAD, 3 = NSGAII, 4 = NSGAIII. Default = SMPSO. Only used if the randomOptimisers parameter is set to 0.
% randomOptimisers --- 1 --- 1 = use random optimisers (SMPSO,MOEAD,NSGAII,NSGAIII) in first phase (defined by delta), and NSGAIII in second phase. 0 = use only specified optimiser. Default = 1 

%------------------------------- Reference --------------------------------
% H. Zille, H. Ishibuchi, S. Mostaghim, and Y. Nojima, A framework for
% large-scale multiobjective optimization based on problem transformation,
% IEEE Transactions on Evolutionary Computation, 2018, 22(2): 260-275.
% ----------------------------------------------------------------------- 
%  Copyright (C) 2020 Heiner Zille
%
%  This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
%  International License. (CC BY-NC-SA 4.0). To view a copy of this license, 
%  visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or see the 
%  pdf-file "License-CC-BY-NC-SA-4.0.pdf" that came with this code. 
%
%  You are free to: 
%  * Share ? copy and redistribute the material in any medium or format
%  * Adapt ? remix, transform, and build upon the material 
%  Under the following terms:
%  * Attribution ? You must give appropriate credit, provide a link to the 
%     license, and indicate if changes were made. You may do so in any reasonable 
%     manner, but not in any way that suggests the licensor endorses you or your use.
%  * NonCommercial ? You may not use the material for commercial purposes.
%  * ShareAlike ? If you remix, transform, or build upon the material, you must 
%    distribute your contributions under the same license as the original.
%  * No additional restrictions ? You may not apply legal terms or technological 
%    measures that legally restrict others from doing anything the license permits.
% 
%  Author of this Code: 
%   Heiner Zille <heiner.zille@ovgu.de> or <heiner.zille@gmail.com>
%
%  This code is based on the following publications:
% 
%  1) Heiner Zille 
%     "Large-scale Multi-objective Optimisation: New Approaches and a Classification of the State-of-the-Art"  
%     PhD Thesis, Otto von Guericke University Magdeburg, 2019 
%     http://dx.doi.org/10.25673/32063 
% 
%  2) Heiner Zille and Sanaz Mostaghim
%     "Comparison Study of Large-scale Optimisation Techniques on the LSMOP Benchmark Functions"  
%     IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Honolulu, Hawaii, November 2017
%     https://ieeexplore.ieee.org/document/8280974 
% 
%  3) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "A Framework for Large-scale Multi-objective Optimization based on Problem Transformation"
%     IEEE Transactions on Evolutionary Computation, Vol. 22, Issue 2, pp. 260-275, April 2018.
%     http://ieeexplore.ieee.org/document/7929324
%  
%  4) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "Weighted Optimization Framework for Large-scale Mullti-objective Optimization"
%     Genetic and Evolutionary Computation Conference (GECCO), ACM, Denver, USA, July 2016
%     http://dl.acm.org/citation.cfm?id=2908979
%
%  This file is intended to work with the PlatEMO framework version 2.5. 
%  Date of publication of this code: 06.04.2020 
%  Last Update of this code: 06.04.2020
%  A newer version of this algorithm may be available. Please contact the author 
%  or see http://www.ci.ovgu.de/Research/Codes.html. 
%
% The files may have been modified in Feb 2021 by the authors of the Platemo framework to work with the Platemo 3.0 release. 
% ----------------------------------------------------------------------- 

    methods
        function main(Algorithm,Problem)
            %% Set the default parameters
            [gamma,groups,psi,t1,t2,q,delta,optimiser,randomOptimisers] = Algorithm.ParameterSet(4,2,3,1000,500,Problem.M+1,0.5,1,1);
            WOF_WeightIndividual.Current(Problem);
            
            % The size of the population of weight-Individuals
            transformedProblemPopulationSize = 10;

            % There are different methods to select the xPrime solutions. 
            % Three methods can be chosen. The first one uses the largest Crowding
            % Distance values from the first non-dominated front. The second one
            % uses a tournament selection on the population based on
            % Pareto-dominance and Crowding Distance. The third option is 
            % introduced in publication (1), see above, based on
            % reference directions for the first m+1 chosen solutions and selects
            % random solutions afterwards. If method 3 is chosen, q is always at
            % least Global.M 
            methodToSelectxPrimeSolutions = 3; 

            diceForType = false;
            if randomOptimisers == 1
                diceForType = true;
            end

            if optimiser == 2 || optimiser == 4 || diceForType == true
                [uniW,Problem.N] = UniformPoint(Problem.N,Problem.M);
            end

            restart = true;


            %% Generate random population
            Population = Problem.Initialization();
            Algorithm.NotTerminated(Population);

            %% Start the alternating optimisation 
            while Problem.FE < delta*Problem.maxFE

                if (diceForType == true)
                    optimiser = randomizeType(optimiser);
                end

                % In case the population is not full (because the number of non-dominated 
                % solutions returned in the previous iteration is smaller than the 
                % population size), it is filled with new solutions. 
                Population = WOFfillPopulation(Population, Problem);

                % Normal optimisation step for t1 evaluations
                if      optimiser == 4
                    Population = WOF_optimiseByNSGAIII(Problem, Population, uniW, t1, false);
                elseif  optimiser == 3
                    Population = WOF_optimiseByNSGAII(Problem, Population, t1, false);
                elseif  optimiser == 2
                    Population = WOF_optimiseByMOEAD(Problem, Population, uniW, t1, false);
                else
                    Population = WOF_optimiseBySMPSO(Problem, Population, t1, false);
                end
                Algorithm.NotTerminated(Population); 

                % Selection of xPrime solutions 
                xPrimeList = WOF_selectxPrimes(Population, q, methodToSelectxPrimeSolutions); 
                WList   = [];

                % do for each xPrime
                for c = 1:size(xPrimeList,2)
                    xPrime              = xPrimeList(c);

                    % create variable groups 
                    [G,gamma]                   = WOF_createGroups(Problem,gamma,xPrime,groups);

                    if (diceForType == true)
                        optimiser = randomizeType(optimiser);
                    end

                    % a dummy object is needed to simulate the global class. Its
                    % necessary to include this method into the Platemo
                    % framework. 
                    GlobalDummy         = WOFcreateGlobalDummy(gamma, xPrime, G, Problem, transformedProblemPopulationSize, psi,optimiser);

                    % Create initial population for the transformed problem
                    WeightPopulation    = WOFcreateInitialWeightPopulation(GlobalDummy.N, gamma, GlobalDummy);

                    % Optimise the transformed problem 
                    if      optimiser == 4
                        WeightPopulation    = WOF_optimiseByNSGAIII(GlobalDummy, WeightPopulation, GlobalDummy.uniW, t2-transformedProblemPopulationSize, true);
                    elseif  optimiser == 3
                        WeightPopulation    = WOF_optimiseByNSGAII(GlobalDummy, WeightPopulation, t2-transformedProblemPopulationSize, true);
                    elseif  optimiser == 2
                        WeightPopulation    = WOF_optimiseByMOEAD(GlobalDummy, WeightPopulation, GlobalDummy.uniW, t2-transformedProblemPopulationSize, true);
                    else 
                        WeightPopulation    = WOF_optimiseBySMPSO(GlobalDummy, WeightPopulation, t2-transformedProblemPopulationSize, true);
                    end

                    % Extract the population 
                    W                   = WOFextractPopulation(WeightPopulation, Problem, Population, G, psi, xPrime, q, methodToSelectxPrimeSolutions);
                    WList               = [WList,W];  
                end

                % Join populations. Duplicate solution (e.g. found in different
                % optimisation steps with different xPrimes) need to be removed. 
                Population          = WOFeliminateDuplicates([Population,WList]);
                Population          = WOFfillPopulation(Population, Problem);

                % Environmental Selection
                [Population,~,~]    = WOF_EnvironmentalSelection(Population,Problem.N);
                Algorithm.NotTerminated(Population);
            end

            %% Optimise until end for uniformity. 
            remainingEvaluations = Problem.maxFE-Problem.FE;

            if ~restart
                % in this case all remaining evaluations are used directly with the
                % secified algorithm
                if (diceForType == true)
                    optimiser = 4;
                end

                if      optimiser == 4
                    Population = WOF_optimiseByNSGAIII(Problem, Population, uniW, remainingEvaluations, false);
                elseif  optimiser == 3
                    Population = WOF_optimiseByNSGAII(Problem, Population, remainingEvaluations, false);
                elseif  optimiser == 2
                    Population = WOF_optimiseByMOEAD(Problem, Population, uniW, remainingEvaluations, false);
                else
                    Population = WOF_optimiseBySMPSO(Problem, Population, remainingEvaluations, false);
                end
            else
                % in this case the value for t1 is used further to do separate
                % chunks of optimisation. 
                while Algorithm.NotTerminated(Population)
                    Population = WOFfillPopulation(Population, Problem);

                    if (diceForType == true)
                        optimiser = 4;
                    end

                    if      optimiser == 4
                        Population = WOF_optimiseByNSGAIII(Problem, Population, uniW, t1, false);
                    elseif  optimiser == 3
                        Population = WOF_optimiseByNSGAII(Problem, Population, t1, false);
                    elseif  optimiser == 2
                        Population = WOF_optimiseByMOEAD(Problem, Population, uniW, t1, false);
                    else
                        Population = WOF_optimiseBySMPSO(Problem, Population, t1, false);
                    end
                end 
            end
        end
    end
end

function GlobalDummy = WOFcreateGlobalDummy(gamma, xPrime, G, Global, populationSize, psi,optimiser)
    % Creates a dummy object. Needed to simulate the global class. Its
    % necessary to include this method into the Platemo
    % framework. 
    GlobalDummy = {};
    GlobalDummy.lower       = zeros(1,gamma);
    GlobalDummy.upper       = ones(1,gamma).*2.0;
    if or(optimiser == 2,optimiser == 4)
        [uniW,GlobalDummy.N]    = UniformPoint(populationSize,Global.M);
        GlobalDummy.uniW        = uniW;
    else
        GlobalDummy.N           = populationSize;
    end
    GlobalDummy.xPrime      = xPrime;
    GlobalDummy.G           = G;
    GlobalDummy.psi         = psi;
    GlobalDummy.xPrimelower = Global.lower;
    GlobalDummy.xPrimeupper = Global.upper;
    GlobalDummy.isDummy     = true;
    GlobalDummy.Global      = Global;
end

function Population = WOFeliminateDuplicates(input)
    % Eliminates duplicates in the population
    [~,ia,~] = unique(input.objs,'rows');
    Population = input(ia);
end

function Population = WOFfillPopulation(input, Global)
    % Fills the population with mutations in case its smaller than Global.N
    Population = input;
    theCurrentPopulationSize = size(input,2);
    if theCurrentPopulationSize < Global.N
        amountToFill    = Global.N-theCurrentPopulationSize;
        FrontNo         = NDSort(input.objs,inf);
        CrowdDis        = CrowdingDistance(input.objs,FrontNo);
        MatingPool      = TournamentSelection(2,amountToFill+1,FrontNo,-CrowdDis);
        Offspring       = OperatorGA(input(MatingPool));
        Population      = [Population,Offspring(1:amountToFill)];
    end
end

function WeightPopulation = WOFcreateInitialWeightPopulation(N, gamma, GlobalDummy)
    %creates an initial population for the transformed problem
    decs = rand(N,gamma).*2.0;
    WeightPopulation = [];
    for i = 1:N
        solution = WOF_WeightIndividual(decs(i,:),GlobalDummy);
        WeightPopulation = [WeightPopulation, solution];
    end
end

function W = WOFextractPopulation(WeightPopulation, Problem, Population, G, psi, xPrime, q, methodToSelectxPrimeSolutions)
    % Extracts a population of individuals for the original problem based
    % on the optimised weights. 
    % First a selection of M+1 Weight-Individuals is selected and apllied
    % to the whole population each. 
    % Second all Weight-Individuals are applied to the chosen xPrime
    % solution, since they are optimised for it. 
    
    % Step 1
    weightIndividualList = WOF_selectxPrimes(WeightPopulation, q, methodToSelectxPrimeSolutions);
    calc = size(weightIndividualList,2)*size(Population,2);

    PopDec1 = ones(calc,Problem.D);
    count = 1;
    for wi = 1:size(weightIndividualList,2)
        weightIndividual = weightIndividualList(wi);
        weightVars = weightIndividual.dec;
        
        for i = 1:size(Population,2)
            individualVars = Population(i).dec;
            
            x = WOF_transformationFunctionMatrixForm(individualVars,weightVars(G),Problem.upper,Problem.lower, psi);

            PopDec1(count,:) = x;
            count = count + 1;
        end
        
    end
    
    W1 = Problem.Evaluation(PopDec1);

    % Step 2
    PopDec2 = [];
    for wi = 1:size(WeightPopulation,2)
        weightIndividual = WeightPopulation(wi);
        weightVars = weightIndividual.dec;
        
            individualVars = xPrime.dec;
            x = 1:Problem.D;
            for j = 1:Problem.D
                x(j) = WOF_transformationFunction(individualVars(j),weightVars(G(j)),Problem.upper(j),Problem.lower(j), psi);   
            end
            PopDec2 = [PopDec2;x]; 
    end
    W2 = Problem.Evaluation(PopDec2);
    
    W = [W1,W2];
end

function optimiser = randomizeType(oldOptimiser)
    
    %randomize with all 4 optimisation algorithms
    optimiser = randi(4);
end
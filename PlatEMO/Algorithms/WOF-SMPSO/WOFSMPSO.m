function WOFSMPSO(Global)
% <algorithm> <W> 
% Weighted optimization framework enhanced SMPSO
% gamma    --- 4    --- Number of groups
% groups   --- 1    --- Grouping method, 1 = linear, 2 = ordered, 3 = random 
% psi      --- 3    --- Transformation function, 1 = Multiplication, 2 = P-Value, 3 = Interval
% t1       --- 1000 --- Number of evaluations for original problem
% t2       --- 500  --- Number of evaluations for transformed problem
% q        ---      --- The number of chosen solutions to do weight optimisation. If no value is specified, the default value is M+1
% delta    --- 0.5  --- The fraction of function evaluations to use for the alternating weight-optimisation phase

%------------------------------- Reference --------------------------------
% H. Zille, H. Ishibuchi, S. Mostaghim, and Y. Nojima, A framework for
% large-scale multiobjective optimization based on problem transformation,
% IEEE Transactions on Evolutionary Computation, 2018, 22(2): 260-275.
%------------------------------- Copyright --------------------------------
%  WOFSMPSO.m 
%  Copyright (C) 2018 Heiner Zille
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%  Author of this Code: 
%   Heiner Zille <heiner.zille@ovgu.de> 
%
%  This file belongs to the following publications:
%
%  1) Heiner Zille and Sanaz Mostaghim
%     "Comparison Study of Large-scale Optimisation Techniques on the LSMOP Benchmark Functions"  
%     IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Honolulu, Hawaii, November 2017
%     https://ieeexplore.ieee.org/document/8280974 
% 
%  2) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "A Framework for Large-scale Multi-objective Optimization based on Problem Transformation"
%     IEEE Transactions on Evolutionary Computation, Vol. 22, Issue 2, pp. 260-275, April 2018.
%     http://ieeexplore.ieee.org/document/7929324
%  
%  3) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "Weighted Optimization Framework for Large-scale Mullti-objective Optimization"
%     Genetic and Evolutionary Computation Conference (GECCO), ACM, Denver, USA, July 2016
%     http://dl.acm.org/citation.cfm?id=2908979
%
%  Date of publication: 12.10.2018 
%  Last Update: 12.10.2018
%--------------------------------------------------------------------------

    %% Set the default parameters
    [gamma,groups,psi,t1,t2,q,delta] = Global.ParameterSet(4,1,3,1000,500,Global.M+1,0.5);
    
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
    
    %% Generate random population
    Population = Global.Initialization();
    Global.NotTermination(Population);
    
    %% Start the alternating optimisation 
    while Global.evaluated < delta*Global.evaluation
        
        % In case the population is not full (because the number of non-dominated 
        % solutions returned in the previous iteration is smaller than the 
        % population size), it is filled with new solutions. 
        Population = fillPopulation(Population, Global);
        
        % Normal optimisation step for t1 evaluations
        Population = WOF_optimiseBySMPSO(Global, Population, t1, false);
        Global.NotTermination(Population); 
        
        % Selection of xPrime solutions 
        xPrimeList = WOF_selectxPrimes(Population, q, methodToSelectxPrimeSolutions); 
        WList   = [];
        
        % do for each xPrime
        for c = 1:size(xPrimeList,2)
            xPrime              = xPrimeList(c);

            % create variable groups 
            G                   = WOF_createGroups(gamma,xPrime,Global.D, groups);
            
            % a dummy object is needed to simulate the global class. Its
            % necessary to include this method into the Platemo
            % framework. 
            GlobalDummy         = createGlobalDummy(gamma, xPrime, G, Global, transformedProblemPopulationSize, psi);
            
            % Create initial population for the transformed problem
            WeightPopulation    = createInitialWeightPopulation(GlobalDummy.N, gamma, GlobalDummy);
           
            % Optimise the transformed problem 
            WeightPopulation    = WOF_optimiseBySMPSO(GlobalDummy, WeightPopulation, t2-transformedProblemPopulationSize, true);
            
            % Extract the population 
            W                   = extractPopulation(WeightPopulation, Global, Population, G, psi, xPrime, q, methodToSelectxPrimeSolutions);
            WList               = [WList,W];  
        end
        
        % Join populations. Duplicate solution (e.g. found in different
        % optimisation steps with different xPrimes) need to be removed. 
        Population          = eliminateDuplicates([Population,WList]);
        Population          = fillPopulation(Population, Global);
        
        % Environmental Selection
        [Population,~,~]    = WOF_EnvironmentalSelection(Population,Global.N);
        Global.NotTermination(Population);
    end
    
    %% Optimise until end for uniformity. 
    remainingEvaluations    = Global.evaluation-Global.evaluated;
    noOfParts               = floor(remainingEvaluations/1000);
    
    for i = 1:noOfParts
        Population = fillPopulation(Population, Global);
        Population = WOF_optimiseBySMPSO(Global, Population, t1, false);
        Global.NotTermination(Population);
    end 
end

function GlobalDummy = createGlobalDummy(gamma, xPrime, G, Global, populationSize, psi)
    % Creates a dummy object. Needed to simulate the global class. Its
    % necessary to include this method into the Platemo
    % framework. 
    GlobalDummy = {};
    GlobalDummy.lower       = zeros(1,gamma);
    GlobalDummy.upper       = ones(1,gamma).*2.0;
    GlobalDummy.N           = populationSize;
    GlobalDummy.xPrime      = xPrime;
    GlobalDummy.G           = G;
    GlobalDummy.psi         = psi;
    GlobalDummy.xPrimelower = Global.lower;
    GlobalDummy.xPrimeupper = Global.upper;
    GlobalDummy.isDummy     = true;
    GlobalDummy.Global      = Global;
end

function Population = eliminateDuplicates(input)
    % Eliminates duplicates in the population
    [~,ia,~] = unique(input.objs,'rows');
    Population = input(ia);
end

function Population = fillPopulation(input, Global)
    % Fills the population with mutations in case its smaller than Global.N
    Population = input;
    theCurrentPopulationSize = size(input,2);
    if theCurrentPopulationSize < Global.N
        amountToFill    = Global.N-theCurrentPopulationSize;
        FrontNo         = NDSort(input.objs,inf);
        CrowdDis        = CrowdingDistance(input.objs,FrontNo);
        MatingPool      = TournamentSelection(2,amountToFill+1,FrontNo,-CrowdDis);
        Offspring       = GA(input(MatingPool));
        Population      = [Population,Offspring(1:amountToFill)];
    end
end

function WeightPopulation = createInitialWeightPopulation(N, gamma, GlobalDummy)
    %creates an initial population for the transformed problem
    decs = rand(N,gamma).*2.0;
    WeightPopulation = [];
    for i = 1:N
        solution = WOF_WeightIndividual(decs(i,:),GlobalDummy);
        WeightPopulation = [WeightPopulation, solution];
    end
end

function W = extractPopulation(WeightPopulation, Global, Population, G, psi, xPrime, q, methodToSelectxPrimeSolutions)
    % Extracts a population of individuals for the original problem based
    % on the optimised weights. 
    % First a selection of M+1 Weight-Individuals is selected and apllied
    % to the whole population each. 
    % Second all Weight-Individuals are applied to the chosen xPrime
    % solution, since they are optimised for it. 
    
    % Step 1
    weightIndividualList = WOF_selectxPrimes(WeightPopulation, q, methodToSelectxPrimeSolutions);
    PopDec1 = [];
    for wi = 1:size(weightIndividualList,2)
        weightIndividual = weightIndividualList(wi);
        weightVars = weightIndividual.dec;
        
        for i = 1:size(Population,2)
            individualVars = Population(i).dec;
            x = 1:Global.D;
            for j = 1:Global.D
                x(j) = WOF_transformationFunction(individualVars(j),weightVars(G(j)),Global.upper(j),Global.lower(j), psi);   
            end
            PopDec1 = [PopDec1;x];
        end
        
    end
    W1 = INDIVIDUAL(PopDec1);

    % Step 2
    PopDec2 = [];
    for wi = 1:size(WeightPopulation,2)
        weightIndividual = WeightPopulation(wi);
        weightVars = weightIndividual.dec;
        
            individualVars = xPrime.dec;
            x = 1:Global.D;
            for j = 1:Global.D
                x(j) = WOF_transformationFunction(individualVars(j),weightVars(G(j)),Global.upper(j),Global.lower(j), psi);   
            end
            PopDec2 = [PopDec2;x]; 
    end
    W2 = INDIVIDUAL(PopDec2);
    
    W = [W1,W2];
end
classdef LCSA < ALGORITHM
% <multi/many> <real/integer> <large/none>
% Linear combination-based search algorithm
% optimiser --- 3 --- The optimisation method. 1 = SMPSO, 2 = NSGA-II, 3 = NSGA-III. Default = NSGA-III

%------------------------------- Reference --------------------------------
% H. Zille, Large-scale Multi-objective Optimisation: New Approaches and a
% Classification of the State-of-the-Art, PhD Thesis, Otto von Guericke
% University Magdeburg, 2019.
%--------------------------------------------------------------------------
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
%     "Linear Search Mechanism for Multi- and Many-Objective Optimisation"
%     10th International Conference on Evolutionary Multi-Criterion Optimization (EMO 2019), 
%        Lecture Notes in Computer Science, vol 11411. 
%        Deb K. et al. (eds), Springer, Cham, East Lansing, Michigan, USA, March 2019  
%     https://doi.org/10.1007/978-3-030-12598-1_32.
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
            optimiser = Algorithm.ParameterSet(3);
            if optimiser == 1
                LCSA_SMPSO(Algorithm,Problem);
            elseif optimiser == 2
                LCSA_NSGAII(Algorithm,Problem);
            else 
                LCSA_NSGAIII(Algorithm,Problem);
            end
        end
        function LCSA_SMPSO(Algorithm,Problem)
            %% Generate random population
            Population       = Problem.Initialization();
            Pbest            = Population;
            [Gbest,CrowdDis] = LCSA_SMPSOUpdateGbest(Population,Problem.N);
            iter = 1;
            xintervall = LCSA_getNumberOfGenerationsForNormalOptimisation(Problem.maxFE, Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Gbest)
                Population       = LCSA_SMPSOOperator(Problem,Population,Pbest,Gbest(TournamentSelection(2,Problem.N,-CrowdDis)));
                [Gbest,CrowdDis] = LCSA_SMPSOUpdateGbest([Gbest,Population],Problem.N);
                Pbest            = LCSA_SMPSOUpdatePbest(Pbest,Population);
                iter = iter + 1;
                if mod(iter,xintervall) == 0
                    Algorithm.NotTerminated(Gbest);
                    noOfXOptGenerations = LCSA_getNumberOfGenerationsForCoefficientOptimisation(Problem.maxFE, Problem.N);
                    decs                = Population.decs;
                    xlower              = ones(1,Problem.N)*-10;	% -10
                    xupper              = ones(1,Problem.N)*10;     % 10
                    newPopSize          = Problem.N;
                    xDecs               = repmat(xlower,newPopSize,1) + rand(newPopSize,Problem.N).*repmat(xupper-xlower,newPopSize,1); %newPopSize solutions with N vars each.
                    xPop                = Problem.Evaluation(transpose(transpose(decs)*transpose(xDecs)),packDecAndVel(xDecs, zeros(newPopSize,Problem.N)));
                    xPbest              = xPop;
                    [xGbest,xCrowdDis]  = LCSA_SMPSOUpdateGbest(xPop,Problem.N);
                    [Gbest,CrowdDis]    = LCSA_SMPSOUpdateGbest([Gbest,xGbest],Problem.N);
                    for temp = 1:noOfXOptGenerations 
                        Algorithm.NotTerminated(Gbest);
                        [xDecs,xVels] = LCSA_CoefficientSMPSOOperator(xPop,xPbest,xGbest(TournamentSelection(2,Problem.N,-xCrowdDis)),xlower,xupper);
                        newSols       = Problem.Evaluation(transpose(transpose(decs)*transpose(xDecs)),packDecAndVel(xDecs, xVels));
                        xPbest        = LCSA_SMPSOUpdatePbest(xPbest,newSols);
                        [xGbest,xCrowdDis] = LCSA_SMPSOUpdateGbest([xPop,newSols],Problem.N);
                        [Gbest,CrowdDis]   = LCSA_SMPSOUpdateGbest([Gbest,xGbest],Problem.N);
                    end 
                end
            end
            function pack = packDecAndVel(xDecs, xVels)
                noOfSolutions = size(xDecs,1);
                T = arrayfun(@(K) struct('xDecs',xDecs(K,:),'xVel',xVels(K,:)), 1:noOfSolutions, 'UniformOutput',0);
                pack = transpose(horzcat(T{:}));
            end
        end
        function LCSA_NSGAII(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = LCSA_NSGAIIEnvironmentalSelection(Population,Problem.N);
            iter = 1;
            xintervall = LCSA_getNumberOfGenerationsForNormalOptimisation(Problem.maxFE, Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,CrowdDis] = LCSA_NSGAIIEnvironmentalSelection([Population,Offspring],Problem.N);
                iter = iter + 1;
                if mod(iter,xintervall) == 0
                    Algorithm.NotTerminated(Population);
                    noOfXOptGenerations = LCSA_getNumberOfGenerationsForCoefficientOptimisation(Problem.maxFE, Problem.N);
                    decs                = Population.decs;
                    xlower              = ones(1,Problem.N)*-10; %-10
                    xupper              = ones(1,Problem.N)*10; %10
                    newPopSize          = Problem.N;
                    xDecs               = repmat(xlower,newPopSize,1) + rand(newPopSize,Problem.N).*repmat(xupper-xlower,newPopSize,1); %newPopSize solutions with N vars each.
                    xPop                = Problem.Evaluation(transpose(transpose(decs)*transpose(xDecs)), xDecs);
                    [~,xFrontNo,xCrowdDis] = LCSA_NSGAIIEnvironmentalSelection(xPop,Problem.N);
                    [Population,~,~] = LCSA_NSGAIIEnvironmentalSelection([Population,xPop],Problem.N);
                    for temp = 1:noOfXOptGenerations 
                        Algorithm.NotTerminated(Population);
                        xMatingPool = TournamentSelection(2,Problem.N,xFrontNo,-xCrowdDis);
                        xDecs       = LCSA_GA(Problem,xPop(xMatingPool).adds,{1,20,1,20},xlower,xupper);
                        newSols     = Problem.Evaluation(transpose(transpose(decs)*transpose(xDecs)), xDecs);
                        [xPop,xFrontNo,xCrowdDis] = LCSA_NSGAIIEnvironmentalSelection([xPop,newSols],Problem.N);
                        [Population,FrontNo,CrowdDis] = LCSA_NSGAIIEnvironmentalSelection([Population,xPop],Problem.N);
                    end 
                end
            end
        end
        function LCSA_NSGAIII(Algorithm,Problem)
            %% Generate the reference points and random population
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
            iter = 1; 
            xintervall = LCSA_getNumberOfGenerationsForNormalOptimisation(Problem.maxFE, Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = LCSA_NSGAIIIEnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
                iter = iter + 1;
                if mod(iter,xintervall) == 0
                    Algorithm.NotTerminated(Population);
                    noOfXOptGenerations = LCSA_getNumberOfGenerationsForCoefficientOptimisation(Problem.maxFE, Problem.N);
                    decs                = Population.decs;
                    xlower              = ones(1,Problem.N)*-10; %-10
                    xupper              = ones(1,Problem.N)*10; %10
                    newPopSize          = Problem.N;
                    xDecs               = repmat(xlower,newPopSize,1) + rand(newPopSize,Problem.N).*repmat(xupper-xlower,newPopSize,1); %newPopSize solutions with N vars each.
                    xPop                = Problem.Evaluation(transpose(transpose(decs)*transpose(xDecs)), xDecs);
                    Zmin                = min([Zmin;xPop(all(xPop.cons<=0,2)).objs],[],1);
                    Population          = LCSA_NSGAIIIEnvironmentalSelection([Population,xPop],Problem.N,Z,Zmin);
                    for temp = 1:noOfXOptGenerations 
                        Algorithm.NotTerminated(Population);
                        xMatingPool = TournamentSelection(2,Problem.N,sum(max(0,xPop.cons),2));
                        xDecs       = LCSA_GA(Problem,xPop(xMatingPool).adds,{1,20,1,20},xlower,xupper);
                        newSols     = Problem.Evaluation(transpose(transpose(decs)*transpose(xDecs)), xDecs);
                        Zmin        = min([Zmin;newSols(all(newSols.cons<=0,2)).objs],[],1);
                        xPop        = LCSA_NSGAIIIEnvironmentalSelection([xPop,newSols],Problem.N,Z,Zmin);
                        Population  = LCSA_NSGAIIIEnvironmentalSelection([Population,xPop],Problem.N,Z,Zmin);
                    end 
                end
            end
        end
    end
end
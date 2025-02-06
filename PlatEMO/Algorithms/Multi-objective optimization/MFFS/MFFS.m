classdef MFFS < ALGORITHM
% <2023> <multi> <binary> 
% Multiform feature selection

%------------------------------- Reference --------------------------------
% R. Jiao, B. Xue, and M. Zhang. Benefiting from single-objective feature 
% selection to multiobjective feature selection: A multiform approach. 
% IEEE Transactions on Cybernetics, 2023, 53(12): 7773-7786.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao

    methods
        function main(Algorithm,Problem)
            % Parameter setting
            N    = round(Problem.N*0.5);
            N1   = round(Problem.N*0.25);
            N2   = Problem.N - N - N1;
            P    = InitializePopulation(Problem);
            alphaSet = [0.01; 0.99]; %Note: We assume the first objective is the selected feature ratio, and the second objective is the classification error rate. 
            PartitionSet = zeros(2, 2);
            [Population, FrontNo, CrowdDis] = EnvironmentalSelectionMOP(P, N);
            [SubPop1, Fitness1] = EnvironmentalSelectionSOP(P, N1, alphaSet(1,:), min(P.objs, [], 1)); % Store solutions with better classification performance 
            [SubPop2, Fitness2] = EnvironmentalSelectionSOP(P, N2, alphaSet(2,:), min(P.objs, [], 1)); % Store solutions with low selected feature ratio         
            BestFitnessSet      = [SubPop1(1).objs; SubPop2(1).objs];
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Offspring generation for multi-objective task
                MatingPool     = TournamentSelection(2, N, FrontNo, -CrowdDis);
                Offspring      = OffspringReproduction(Problem, Population(MatingPool), [Population, SubPop1, SubPop2]);
                % Offspring generation for single-objective task I
                MatingPoolPop1 = TournamentSelection(2, N1, Fitness1);
                OffPop1        = OffspringReproduction(Problem, SubPop1(MatingPoolPop1), [Population, SubPop1, SubPop2, Offspring]);
                % Offspring generation for single-objective task II
                MatingPoolPop2 = TournamentSelection(2, N2, Fitness2);
                OffPop2        = OffspringReproduction(Problem, SubPop2(MatingPoolPop2), [Population, SubPop1, SubPop2, Offspring, OffPop1]);
                % Environmental selection for multi-objective task
                [Population, FrontNo, CrowdDis] = EnvironmentalSelectionMOP([Population, Offspring, OffPop1, OffPop2], N);
                % Update ideal vector
                Fmin = min(Population.objs, [], 1);
                % Environmental selection for single-objective task I
                [SubPop1, Fitness1] = EnvironmentalSelectionSOP([Offspring, OffPop2, SubPop1, OffPop1], N1, alphaSet(1,:), Fmin); 
                % Environmental selection for single-objective task II
                [SubPop2, Fitness2] = EnvironmentalSelectionSOP([Offspring, OffPop1, SubPop2, OffPop2], N2, alphaSet(2,:), Fmin);                
                % Judge whether single-objective tasks convergenced
                [flag, BestFitnessSet] = boolImprovement(BestFitnessSet, [SubPop1(1).objs; SubPop2(1).objs], 5);
                % Calculate direction vectors
                [alphaSet, PartitionSet] = CalWeight(PartitionSet, alphaSet, Population, FrontNo, flag, Fmin);
                if flag==1
                    [SubPop1, Fitness1, SubPop2, Fitness2] = ReInitialization(Problem, Population, N1, N2, alphaSet, Fmin);
                end
                % Population = FSTestCOP(Problem, Population); % Apply to the test set
            end
        end
    end
end


function [alphaSet, PartitionSet] = CalWeight(PartitionSet, alphaSet, Population, FrontNo, flag, Fmin)
    % Calculate direction vectors based on partion point and ideal point
    if flag == 1
        PartitionSet = CalPartitionPoint(Population, FrontNo, PartitionSet);
        alphaSet = UpdateWeight(PartitionSet, Fmin);
    else
        if ~isequal(alphaSet, [0.01; 0.99])
            alphaSet = UpdateWeight(PartitionSet, Fmin);
        end
    end
end

function alphaSet = UpdateWeight(PartitionSet, Fmin)
    % Update direction vectors
    if isequal(PartitionSet(1, :), PartitionSet(2, :)) && ~isequal(PartitionSet(1, :), zeros(1, 2))
        Alpha1 = round(rand(), -2);
        Alpha2 = (PartitionSet(2, 1) - Fmin(1, 1))./((PartitionSet(2, 1) - Fmin(1, 1)) + (PartitionSet(1, 2) - Fmin(1, 2)));
        Alpha2 = round(Alpha2, -2);
    elseif isequal(PartitionSet(1, :), zeros(1, 2)) && isequal(PartitionSet(2, :), zeros(1, 2))
        Alpha1 = round(rand(), -2);
        Alpha2 = round(rand(), -2);
    else
        Alpha1 = (PartitionSet(1, 1)-Fmin(1,1))./((PartitionSet(1, 1) - Fmin(1, 1)) + (PartitionSet(1, 2) - Fmin(1, 2)));
        Alpha1 = round(Alpha1, -2);
        Alpha2 = (PartitionSet(2, 1) - Fmin(1, 1))./((PartitionSet(2, 1) - Fmin(1, 1)) + (PartitionSet(2, 2) - Fmin(1, 2)));
        Alpha2 = round(Alpha2, -2);
    end
    alphaSet = [Alpha1; Alpha2];
end


function Population = InitializePopulation(Problem)
    % Initialization for all tasks
    Pop = zeros(Problem.N, Problem.D);
    for i = 1 : Problem.N
        k = randperm(round(Problem.D), 1);
        j = randperm(Problem.D, k);
        Pop(i, j) = 1;
    end
    Population = Problem.Evaluation(Pop);
end

function [SubPop1, Fitness1, SubPop2, Fitness2] = ReInitialization(Problem, Pop, N1, N2, alphaSet, Fmin)
    % Reinitialize subpopulations for the single-objective task
    Pop = Pop.decs;
    for i =1:size(Pop, 1)
            index1 = find( Pop(i, :));
            index2 = find(~Pop(i, :));
            if size(index1, 2) > 0
                Pop(i, index1(randi(end, 1, 1))) = 0;
            end
            if size(index2, 2) > 0
                Pop(i, index2(randi(end, 1, 1))) = 1;
            end
    end
    Population = Problem.Evaluation(Pop);
    [SubPop1, Fitness1] = EnvironmentalSelectionSOP(Population, N1, alphaSet(1,:), Fmin); 
    [SubPop2, Fitness2] = EnvironmentalSelectionSOP(Population, N2, alphaSet(2,:), Fmin);
end

function [flag, BestFitnessSet] = boolImprovement(BestFitnessSet, bestF, gen)
    % Judge whether the best fitness has improved in gen generations
    flag = 0;
    if size(BestFitnessSet, 2)/2 < gen
        BestFitnessSet = [BestFitnessSet, bestF];
    else
        for j=1:gen - 1
            BestFitnessSet(:, 2*j-1:2*j) = BestFitnessSet(:, 2*(j+1)-1:2*(j+1));
        end
        BestFitnessSet(:, 2*gen-1:2*gen) = bestF;
        if isequal(BestFitnessSet(:,1:2), BestFitnessSet(:,2*gen-1:2*gen))
            flag = 1;
            BestFitnessSet = [];
        end
    end
end

function fitness = FitnessSOP(Population, alpha, z)
    % The fitness function of single-objective task
    Obj     = Population.objs;
    Obj1    = Obj(:, 1); % Selected feature ratio
    Obj2    = Obj(:, 2); % Classification error rate
    fitness = (Obj1 - z(:, 1)).*(1-alpha) + (Obj2 - z(:, 2)).*alpha;
end

function [Population, fitness] = EnvironmentalSelectionSOP(Population, N, alpha, z)
    % Environmental selection for single-objective task
    Fit        = FitnessSOP(Population, alpha, z);
    [~, Rank]  = sort(Fit, 'ascend');
    Population = Population(Rank(1:N));
    fitness    = Fit(Rank(1:N));
end

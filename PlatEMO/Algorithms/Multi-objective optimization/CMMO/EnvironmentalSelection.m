function [Population,Fitness,D_Dec,D_Pop] = EnvironmentalSelection(Population,N)
% The environmental selection of SPEA2 based on objective space

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    %% Calculate the fitness of each solution
    [D_Dec,D_Pop,Fitness] = CalFitness(Population.objs,Population.decs);

    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).objs,Population(Next).decs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
    D_Dec      = D_Dec(Next);
    D_Pop      = D_Pop(Next);
    % Sort the population
    [Fitness,rank] = sort(Fitness);
    Population = Population(rank);
    D_Dec = D_Dec(rank);
    D_Pop = D_Pop(rank);
end

function Del = Truncation(PopObj,PopDec,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance_Pop = pdist2(PopObj,PopObj);
    Distance_Pop(logical(eye(length(Distance_Pop)))) = inf;
    
    Distance_Dec = pdist2(PopDec,PopDec);
    Distance_Dec(logical(eye(length(Distance_Dec)))) = inf;
    
    D = Distance_Pop + Distance_Dec;
    
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(D(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
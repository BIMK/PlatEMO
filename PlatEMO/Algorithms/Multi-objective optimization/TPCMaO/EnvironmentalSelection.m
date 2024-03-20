function [Population,Fitness,feasible_rate] = EnvironmentalSelection(Population,N,type,epsilon)
% Three environmental selection strategies

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the fitness of each solution
    switch type
        case 1
            Fitness = CalFitness(Population.objs,Population.cons);
        case 2
            Fitness = CalFitness(Population.objs);
        case 3
            PopCon = Population.cons;
            CV     = sum(max(0,PopCon),2);
            index  = find(CV<=epsilon);
            PopCon(index(:,1),:) = 0;
            Fitness = CalFitness(Population.objs,PopCon);
    end
    [~,extreme] = min(Population.objs,[],1);
    Fitness(unique(extreme)) = 0;
    
    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).objs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
    CV = sum(max(0,Population.cons),2);
    feasible_rate = sum(CV==0)/size(Population,2);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    N = size(PopObj,1);
    
    %% Calculate the shifted distance between each two solutions
    Distance = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    
    %% Truncation
    Del = false(1,N);
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
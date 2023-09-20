function [Population,Fitness,selected] = EnvironmentalSelection(Population,N,processcon,totalcon,epsilon)
% The environmental selection of SPEA2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the fitness of each solution
    if nargin == 5
        Fitness = CalFitness1(Population.objs,Population.cons,processcon,epsilon);
    else
        if ~processcon
            Fitness = CalFitness(Population.objs);
        elseif processcon > totalcon
            Fitness = CalFitness(Population.objs,Population.cons);
        else
            Fitness = CalFitness(Population.objs,Population.cons,processcon);
        end
    end
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
    selected(1) = size(find(Next(end-N+1:end) == 1),2); %上一代选中
    selected(2) = size(find(Next(1:N) == 1),2); %上一代生成的后代选中
    p = (size(Next,2) - 2*N ) / (N/2);
    for j = 1 : p
        selected(j+2) = size(find(Next(N+1+((j-1)*N/2):N+j*N/2) == 1),2);
    end
    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
    % Sort the population
    [Fitness,rank] = sort(Fitness);
    Population = Population(rank);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
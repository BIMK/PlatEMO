function Population = UpdatePopulation(Population,New,N,status)
% Select NI-mu solutions from Population and combine them with mu new solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete the solutions in Population that may duplicate with new solutions
    [~,index]  = setdiff(Population.objs,New.objs,'rows');
    Population = Population(index);
    
    %% Delete the duplicated solutions in current Population
    [~, Unduplicated] = unique(Population.objs,'rows');
    Population        = Population(Unduplicated);

    %% Calculate the fitness of each solution  
    if nargin == 3
        Fitness = CalFitness(Population.objs);
    else
        if status ==1
            Fitness = CalFitness(Population.objs,Population.cons);
        elseif status ==2
            CV      = sum(max(0,Population.cons),2);         
            Fitness = CalFitness([Population.objs,CV]);
        elseif status ==3
            Fitness = CalFitness(Population.objs,Population.cons);
        end
    end
    
    %% Environmental selection 
    if length(Population)>N
        Next = Fitness < 1;
        if sum(Next) < N
            [~,Rank] = sort(Fitness);
            Next(Rank(1:N)) = true;
        elseif sum(Next) > N
            Del  = Truncation(Population(Next).objs,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
        Population = Population(Next);
    end
    Population = cat(2,Population,New);  
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
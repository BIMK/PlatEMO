function Population = subSPEA2(Population,N)
% The environmental selection of SPEA2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    %% Calculate the fitness of each solution
    Fitness = CalFitness(Population.objs);
    CC = min([length(Population),N]);
    
    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < CC
        [~,Rank] = sort(Fitness);
        Next(Rank(1:CC)) = true;
    elseif sum(Next) > CC
        Del  = Truncation(Population(Next).objs,sum(Next)-CC);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
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

function Fitness = CalFitness(PopObj)
% Calculate the fitness of each solution
    N = size(PopObj,1);

    %% Detect the dominance relation between each two solutions
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end
    
    %% Calculate S(i)
    S = sum(Dominate,2);
    
    %% Calculate R(i)
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end
    
    %% Calculate D(i)
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D = 1./(Distance(:,floor(sqrt(N)))+2);
    
    %% Calculate the fitnesses
    Fitness = R + D';
end
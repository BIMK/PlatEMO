function Population = Archive(Population,N)
% Select feasible and non-dominated solutions by using SPEA2-CDP

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    fIndex     = all(Population.cons <= 0,2);
    Population = Population(fIndex);
    if isempty(Population)
        return;
    elseif size(Population,2) > N
        Fitness = CalFitness(Population.objs,Population.cons);
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
    end
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
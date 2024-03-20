function Population = archive(Population,N)
% Select feasible and non-dominated solutions by SPEA2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select feasible solutions
    fIndex           = all(Population.cons <= 0,2);
    Population       = Population(fIndex);
    if isempty(Population)
        return
    else
        Fitness = CalFitness(Population.objs);
        Next    = Fitness < 1;
        if sum(Next) > N
            Del  = Truncation(Population(Next).objs,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
        Population = Population(Next);
    end
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
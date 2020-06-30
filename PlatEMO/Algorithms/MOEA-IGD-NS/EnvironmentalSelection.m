function Population = EnvironmentalSelection(Population,W,N)
% The environmental selection of MOEA/IGD-NS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;
    
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    Choose = LastSelection(Population(Last).objs,W,N-sum(Next));
    Next(Last(Choose)) = true;
    % Population for next generation
    Population = Population(Next);
end

function Remain = LastSelection(PopObj,W,K)
% Select part of the solutions in the last front

    N  = size(PopObj,1);
    NW = size(W,1);

    %% Calculate the distance between each solution and point
    Distance = pdist2(PopObj,W);
    Con      = min(Distance,[],2);
    
    %% Delete the solution which has the smallest metric contribution one by one
    [dis,rank] = sort(Distance,1);
    Remain     = true(1,N);
    while sum(Remain) > K
        % Calculate the fitness of outliers
        Outliers = Remain;
        Outliers(rank(1,:)) = false;
        METRIC   = sum(dis(1,:)) + sum(Con(Outliers));
        Metrics  = inf(1,N);
        Metrics(Outliers) = METRIC - Con(Outliers);
        % Calculate the fitness of other solutions
        for p = find(Remain & ~Outliers)
            temp = rank(1,:) == p;
            outliers = false(1,N);
            outliers(rank(2,temp)) = true;
            outliers = outliers & Outliers;
            Metrics(p) = METRIC - sum(dis(1,temp)) + sum(dis(2,temp)) - sum(Con(outliers));
        end
        % Delete the worst solution and update the variables
        [~,del] = min(Metrics);
        temp = rank ~= del;
        dis  = reshape(dis(temp),sum(Remain)-1,NW);
        rank = reshape(rank(temp),sum(Remain)-1,NW);
        Remain(del) = false;
    end
end
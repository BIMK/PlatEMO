function [Population,Range] = EnvironmentalSelection(Population,RefPoint,Range,N)
% The environmental selection of AR-MOEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    CV = sum(max(0,Population.cons),2);
    if sum(CV==0) > N
        %% Selection among feasible solutions
        Population = Population(CV==0);
        % Non-dominated sorting
        [FrontNo,MaxFNo] = NDSort(Population.objs,N);
        Next = FrontNo < MaxFNo;
        % Select the solutions in the last front
        Last   = find(FrontNo==MaxFNo);
        Choose = LastSelection(Population(Last).objs,RefPoint,Range,N-sum(Next));
        Next(Last(Choose)) = true;
        Population = Population(Next);
        % Update the range for normalization
        Range(2,:) = max(Population.objs,[],1);
        Range(2,Range(2,:)-Range(1,:)<1e-6) = 1;
    else
        %% Selection including infeasible solutions
        [~,rank]   = sort(CV);
        Population = Population(rank(1:N));
    end
end

function Remain = LastSelection(PopObj,RefPoint,Range,K)
% Select part of the solutions in the last front

    N  = size(PopObj,1);
    NR = size(RefPoint,1);

    %% Calculate the distance between each solution and point
    Distance    = CalDistance(PopObj-repmat(Range(1,:),N,1),RefPoint);
    Convergence = min(Distance,[],2);
    
    %% Delete the solution which has the smallest metric contribution one by one
    [dis,rank] = sort(Distance,1);
    Remain     = true(1,N);
    while sum(Remain) > K
        % Calculate the fitness of noncontributing solutions
        Noncontributing = Remain;
        Noncontributing(rank(1,:)) = false;
        METRIC = sum(dis(1,:)) + sum(Convergence(Noncontributing));
        Metric = inf(1,N);
        Metric(Noncontributing) = METRIC - Convergence(Noncontributing);
        % Calculate the fitness of contributing solutions
        for p = find(Remain & ~Noncontributing)
            temp = rank(1,:) == p;
            noncontributing = false(1,N);
            noncontributing(rank(2,temp)) = true;
            noncontributing = noncontributing & Noncontributing;
            Metric(p) = METRIC - sum(dis(1,temp)) + sum(dis(2,temp)) - sum(Convergence(noncontributing));
        end
        % Delete the worst solution and update the variables
        [~,del] = min(Metric);
        temp    = rank ~= del;
        dis     = reshape(dis(temp),sum(Remain)-1,NR);
        rank    = reshape(rank(temp),sum(Remain)-1,NR);
        Remain(del) = false;
    end
end
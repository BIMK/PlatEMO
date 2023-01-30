function [Population,FrontNo] = EnvironmentalSelection(Population,N)
% The environmental selection of OSP-NSDE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Elaine Guerrero-Pena

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;

    %% Select the solutions in the last front by truncation
    if sum(Next) < N
        Last = find(FrontNo==MaxFNo);
        apt = Population.objs;
        Pop = apt(Last,:);
        Next(Last) = true;
        Temp = Last;
    elseif sum(Next) > N
        apt = Population.objs;
        Pop = apt(Next,:);
        Temp = find(Next);
    end
    Del  = Truncation(Pop,sum(Next)-N);
    Next(Temp(Del)) = false;

    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
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

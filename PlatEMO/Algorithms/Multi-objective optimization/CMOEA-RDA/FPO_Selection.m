function FPO = FPO_Selection(FPO, Offspring, N)
%The environment selection method of FPO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    %% Dominance relation calculation
    FPO          = [FPO, Offspring];
    NonDominated = DominationCal(FPO, 1);
    FPO          = FPO(NonDominated);
    PopObj       = FPO.objs;
    [P, ~]       = size(PopObj);
    if P <= N
        FPO = FPO([1:P, unidrnd(P,[1, N - P])]);
    else
        PopCon = FPO.cons;
        Cons   = sum(max(0,PopCon),2);
        if sum(Cons <= 0) <= N
            [~, index]        = sort(Cons);
            remained_solution = index(1: N);
            FPO               = FPO(remained_solution);
        else
            FPO    = FPO(Cons <= 0);
            PopObj = FPO.objs;
            [P, ~] = size(PopObj);
            zmin   = min(PopObj, [], 1) - 1e-6;
            zmax   = max(PopObj, [], 1);
            PopObj = (PopObj - zmin) ./ (zmax - zmin);
            Del    = Truncation(PopObj, P - N);
            remained_solution      = 1 : P;
            remained_solution(Del) = [];
            FPO = FPO(remained_solution);
        end
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
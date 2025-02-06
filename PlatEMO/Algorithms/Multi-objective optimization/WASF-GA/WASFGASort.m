function [FrontNo,MaxFNo] = WASFGASort(Vectors, PopObj, nsort, Point, ro)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [nvectors,~] = size(Vectors);
    [Loc,MaxFNo] = frontsclass(Vectors, PopObj,inf, Point, ro);
    [popsize,~]  = size(PopObj);
    FrontNo      = inf(1,size(PopObj,1));

    for i = 1 : popsize
        Position = find(Loc == i);
        iter = 0;
        while nvectors*iter < Position
            iter = iter + 1;
        end
        if iter == 0 || iter > nsort
            FrontNo(i) = inf;
        else
            FrontNo(i) = iter;
        end
    end
end

function [Loc, Max] = frontsclass(Vectors, PopObj, nsort, Point, ro)
    [nvectors,~] = size(Vectors);
    [N,~]  = size(PopObj);
    FrontG = [];
    % SolG will store the different solutions sorted by the achievement
    % scalarizing function
    SolG    = [];
    PopObj2 = PopObj;
    Max     = 0;
    while length(FrontG) < N && Max < nsort
        % n will be the size of the population that will be compared in
        % each iteration, it will change in every iteration.
        [n,~] = size(PopObj);
        Max   = Max + 1;
        for i = 1 : min(nvectors, n)
            Front  = [];
            Values = zeros(n, 1);
            for j = 1 : n
                Values(j) = max((PopObj(j, :) - Point) .* Vectors(i, :)) + ro * sum(Vectors(i, :) .* (PopObj(j, :) - Point));
            end
            Vmin      = min(Values);
            Sol1      = find(Values == Vmin);
            Front     = [Front, Sol1(1)];
            FrontG    = [FrontG, Front];
            valid_Loc = Front(Front <= size(PopObj, 1));
            SolG      = [SolG; PopObj(valid_Loc, :)];
            PopObj(valid_Loc,:) = [];
            [n,~]     = size(PopObj);
        end
    end
    Loc = find_Loc(SolG, PopObj2);
end

function index = find_Loc(rows, initial_matrix)
% Identify the position in the original matrix of the solutions in SolG

    index = zeros(size(rows, 1), 1);
    for i = 1 : size(rows, 1)
        [equalrow,loc] = ismember(rows(i, :), initial_matrix, 'rows');
        if equalrow
            index(i) = loc;
        end
    end
end
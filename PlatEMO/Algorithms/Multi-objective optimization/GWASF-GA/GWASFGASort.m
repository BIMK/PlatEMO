function [FrontNo,MaxFNo] = GWASFGASort(Vectors, PopObj, Utop,Nadir, nsort, ro, eps)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [nvectors,~] = size(Vectors);
    [Loc,MaxFNo] = frontloc(Vectors, PopObj,Utop, Nadir, inf, ro, eps);
    [popsize,~]  = size(PopObj);
    FrontNo = inf(1,size(PopObj,1));
    for i = 1 : popsize
        position = find(Loc == i);
        count = 0;
        while nvectors*count < position
            count = count + 1;
        end
        if count == 0 || count > nsort
            FrontNo(i) = inf;
        else
            FrontNo(i) = count;
        end
    end
end

function [Loc, Max] = frontloc(Vectors, PopObj,Utop,nadir, nsort, ro, eps)
    [lengthVectors,~] = size(Vectors);
    [bound, ~]        = size(PopObj);
    % SolG will store the different solutions sorted by the achievement
    % scalarizing function
    SolutionsG = [];
    PopObj2    = PopObj;
    Max        = 0;
    while size(SolutionsG,1) < bound && Max < nsort
        % n will be the size of the population that will be compared in
        % each iteration, it will change in every iteration.
        [n,~] = size(PopObj);
        Max   = Max + 1;
        % At each iteration, we alternate between the nadir and utopian points to sort the population into frontiers
        for i = 1 : (lengthVectors/2)
            ValuesU = zeros(n, 1);
            for j = 1 : n
              ValuesU(j) = max((PopObj(j, :) - Utop) .* Vectors((2*i -1), :)) + ro * sum(Vectors((2*i -1 ), :) .* (PopObj(j, :) - Utop));
            end
            [~,indexU] = sort(ValuesU);
            Sol1       = indexU(1);
            SolutionsG = [SolutionsG; PopObj(Sol1, :)];
            if size(SolutionsG,1) == bound
                break;
            end
            PopObj(Sol1,:) = [];
            [n,~]   = size(PopObj);
            ValuesN = zeros(n,1);
            for j = 1 : n
              ValuesN(j) = max((PopObj(j, :) - nadir) .* Vectors((2*i ), :)) + ro * sum(Vectors((2*i ), :) .* (PopObj(j, :) - nadir));
            end
            [~,indexN] = sort(ValuesN);
            Sol2       = indexN(1);                              
            SolutionsG = [SolutionsG; PopObj(Sol2, :)];
            if size(SolutionsG,1) == bound
                break;
            end
            PopObj(Sol2,:) = [];
            [n,~] = size(PopObj);
        end
    end
    Loc = find_Loc(SolutionsG, PopObj2);
end

function location = find_Loc(moved_rows, initial_matrix)
% Identify the position in the original matrix of the solutions in SolG

    location = zeros(size(moved_rows, 1), 1);
    for i = 1 : size(moved_rows, 1)
        % Find the position of the moved row in the initial matrix
        [equal_row,index] = ismember(moved_rows(i, :), initial_matrix, 'rows');
        % Verify if there is any similarity
        if equal_row
            location(i) = index;
        end
    end
end
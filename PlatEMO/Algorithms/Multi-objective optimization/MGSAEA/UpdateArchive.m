function Archive = UpdateArchive(Archive,N)
% Update the archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete the duplicated solutions
    [~, Unduplicated] = unique(Archive.objs,'rows');
    Archive           = Archive(Unduplicated);
    if length(Archive) > N
        % Calculate the fitness of each solution
        Fitness = CalFitness(Archive.objs, Archive.cons);          
        % Environmental selection
        Next = Fitness < 1;
        if sum(Next) < N
            [~,Rank] = sort(Fitness);
            Next(Rank(1:N)) = true;
        elseif sum(Next) > N
            Del  = Truncation(Archive(Next).objs,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
        Archive = Archive(Next);
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
function [Archive,sucessful] = UpdateArchive(Archive,Offspring)
% Update the archive by NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = length(Archive);
    Archive = [Archive,Offspring];

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Archive.objs,N);
    Next = FrontNo < MaxFNo;
    
    %% Select the solutions in the last front
    Last = find(FrontNo==MaxFNo);
    CrowdDis = CrowdingDistance(Archive(Last).objs);
    [~,rank] = sort(CrowdDis,'descend');
    Next(Last(rank(1:N-sum(Next)))) = true;
    
    %% Update the archive
    Archive   = Archive(Next);
    sucessful = Next(N+1:end);
end

function CrowdDis = CrowdingDistance(PopObj)
% Calculate the crowding distance of each solution in the same front

    [N,M]    = size(PopObj);
    
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);
    for i = 1 : M
        [~,rank] = sortrows(PopObj(:,i));
        CrowdDis(rank(1))   = inf;
        CrowdDis(rank(end)) = inf;
        for j = 2 : N-1
            CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end
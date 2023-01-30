function [Archive,sucessful] = UpdateArchive(Archive,Offspring)
% Update the archive by NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
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
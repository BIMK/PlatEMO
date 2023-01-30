function Del = EnvironmentalSelection(Population,nSub)
% The environmental selection of IM-MOEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,nSub);
    Next = FrontNo < MaxFNo;

    %% Select the solutions in the last front by crowding distance
    Last     = find(FrontNo==MaxFNo);
    CrowdDis = CrowdingDistance(Population(Last).objs);
    [~,rank] = sort(CrowdDis,'descend');
    Next(Last(rank(1:nSub-sum(Next)))) = true;
    % The index of deleted solutions
    Del = ~Next;
end
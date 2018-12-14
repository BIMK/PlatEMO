function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N,Psi)
% The environmental selection of NSGA-II-conflict

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Selected = zeros(1,N);
    FrontNo  = zeros(1,N);
    CrowdDis = zeros(1,N);
    PopObj   = Population.objs;
    for i = 1 : length(Psi)
        index = (i-1)*ceil(N/length(Psi))+1 : min(N,i*ceil(N/length(Psi)));
        [Selected(index),FrontNo(index),CrowdDis(index)] = SubSelection(PopObj(:,Psi{i}),length(index));
    end
    Population = Population(Selected);
end

function [Next,FrontNo,CrowdDis] = SubSelection(PopObj,N)
% Environmental selection based on only several objectives

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(PopObj,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(PopObj,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Next     = find(Next);
    FrontNo  = FrontNo(Next);
    CrowdDis = CrowdDis(Next);
end
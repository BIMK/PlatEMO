function [Population,FrontNo,CrowdDis] = DominationSelection(Global,Population)
% The dominant relationship and crowding based environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Shufen Qin
% E-mail: shufen.qin@stu.tyust.edu.cn

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Global.N);
    Next = false(1,length(FrontNo));
    Next(FrontNo<MaxFNo) = true;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:Global.N-sum(Next)))) = true;
    
    %% Population for next generation
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
    Population = Population(Next);
end
function [PopDec,PopObj,FrontNo,CrowdDis] = EnvironmentalSelection(PopDec,PopObj,M,N)
% The environmental selection of M1-2 based on non-dominating front number and
% crowded distance with the assistance of surrogate models

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete the duplicated points 
    [~, Unduplicated] = unique(PopObj(:,1:M),'rows');
    PopDec            = PopDec(Unduplicated,:);
    PopObj            = PopObj(Unduplicated,:);
    
    %% Non-dominated sorting
    RealObj = PopObj(:,1:M);
    PopCon  = PopObj(:,M+1:end);
    CV      = sum(max(0,PopCon),2);
    [FrontNo,MaxFNo] = NDSort(RealObj,CV,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(RealObj,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    if length(Rank) >= N-sum(Next)
        Next(Last(Rank(1:N-sum(Next)))) = true;
    else
        Next(Last(Rank(1:length(Rank)))) = true;
    end
    
    %% Population for next generation
    PopDec     = PopDec(Next,:);
    PopObj     = PopObj(Next,:);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);    
end
function [Population,FrontNo,d2] = EnvironmentalSelection(OffSpring,W,N)
% The selection of TSNSGAII

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    %% Normalization
    PopObj = OffSpring.objs;
    Fmin   = min(PopObj,[],1);
    Fmax   = max(PopObj,[],1);
    PopObj = (PopObj-repmat(Fmin,size(PopObj,1),1))./repmat(Fmax-Fmin,size(PopObj,1),1);
    
    %% Association
    normP   = sqrt(sum(PopObj.^2,2));
    Cosine  = 1 - pdist2(PopObj,W,'cosine');
    d1      = repmat(normP,1,size(W,1)).*Cosine;
    d2      = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2);
    [d2,RP] = min(d2,[],2);
    d1      = d1((1:length(RP))'+(RP-1)*length(RP));
    
    %% Favor extreme solutions
    ND              = find(NDSort(PopObj,1)==1);
    [~,Extreme]     = max(PopObj(ND,:),[],1);
    d1(ND(Extreme)) = 0;
    d2(ND(Extreme)) = 0;
    
    %% SPDsorting
    [FrontNo,MaxFNo] = SPDSort(PopObj,d1,d2,RP,N);
    Next = FrontNo < MaxFNo;

    %% Select the solutions in the last front
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(d2(Last));
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = OffSpring(Next);
    FrontNo    = FrontNo(Next);
    d2         = d2(Next);
end
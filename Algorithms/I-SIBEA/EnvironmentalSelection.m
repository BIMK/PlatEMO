function [Population,FrontNo,WHVLoss] = EnvironmentalSelection(Population,N,wz,AA,RA)
% The environmental selection of I-SIBEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the WHV loss of each solution front by front
    WHVLoss = CalWHVLoss(Population.objs,FrontNo,wz,AA,RA);
    
    %% Select the solutions in the last front based on their WHV loss
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(WHVLoss(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    WHVLoss    = WHVLoss(Next);
end
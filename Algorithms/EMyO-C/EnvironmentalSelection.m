function [Population] = EnvironmentalSelection(Population,N,Z)
% The environmental selection of EMyOC

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    
    Next = FrontNo<MaxFNo;
    Last = find(FrontNo==MaxFNo);
    
    %% Select the solutions in the last front using clustering based truncation
    index = Truncation(Population(Last).objs,Z,N-sum(Next));  
    
    Next(Last(index)) = true;
    
    %% Population for next generation
    Population = Population(Next);
end
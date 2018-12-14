function Population = EnvironmentalSelection(Population,N,objective_set)
% The environmental selection of SIBEA-kEMOSS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Liangli Zhen

    PopObj = Population.objs;
    
    %% Non-dominated sorting    
    [FrontNo,MaxFNo] = NDSort(PopObj(:,objective_set),N);
    Next = FrontNo < MaxFNo;
    
    %% Select the solutions in the last front based on their HV loss
    Last = find(FrontNo==MaxFNo);
    % Calculate the WHV loss of each solution in the last front
    HVLoss   = CalHVLoss(PopObj(Last,objective_set),FrontNo(Last));
    [~,Rank] = sort(HVLoss,'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
end
function Population = EnvironmentalSelection(Population,N)
% The environmental selection of BiGE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting wrt the actual objectives
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;
    
    %% Proximity and crowding degree estimation for the last-front solutions of the actual objectives
    Last  = find(FrontNo==MaxFNo);
    BiObj = Estimation(Population(Last).objs,1/N^(1/length(Population(1).obj)));
    
    %% Non-dominated sorting wrt the bi-criteria
    [FrontNo2,MaxFNo2] = NDSort(BiObj,N-sum(Next));
    Next(Last(FrontNo2<MaxFNo2)) = true;

    %% Select the solutions in the last front
    Last2 = Last(FrontNo2==MaxFNo2);
    Next(Last2(randperm(length(Last2),N-sum(Next)))) = true;
    % Population for next generation
    Population = Population(Next);
end
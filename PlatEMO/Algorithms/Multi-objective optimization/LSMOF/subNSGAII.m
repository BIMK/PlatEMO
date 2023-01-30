function Population = subNSGAII(Problem,Population,Operator,N)
% Sub-optimizer in LSMOF (NSGA-II)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He
    
    FrontNo  = NDSort(Population.objs,inf);
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    if Operator == 1
        MatingPool = TournamentSelection(2,N,FrontNo,-CrowdDis);
        Offspring  = OperatorGA(Problem,Population(MatingPool));
    else
        MatingPool1 = TournamentSelection(2,N,FrontNo,-CrowdDis);
        MatingPool2 = TournamentSelection(2,N,FrontNo,-CrowdDis);
        MatingPool3 = TournamentSelection(2,N,FrontNo,-CrowdDis);
        Offspring   = OperatorDE(Problem,Population(MatingPool1),Population(MatingPool2),Population(MatingPool3));
    end
    Population = EnvironmentalSelection([Population,Offspring],N);
end
function MatingPool = MatingSelection(PopObj)
% The mating selection of AGE-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the front number and crowding distance of each solution
    FrontNo  = NDSort(PopObj,inf);
    CrowdDis = CrowdingDistance(PopObj,FrontNo);
    
    %% Reduce the population
    Remain = find(rand(1,size(PopObj,1))<1./FrontNo);

    %% Binary tournament selection
    MatingPool = Remain(TournamentSelection(2,size(PopObj,1),-CrowdDis(Remain)));
end
function [Pgood,Pbad] = NDS(Population,K)
% Sort the population based on non-dominated sorting and crowding distance

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    FrontNo  = NDSort(Population.objs,inf);
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    [~,rank] = sortrows([FrontNo;-CrowdDis]');
    Pgood    = Population(rank(1:K));
    Pbad     = Population(rank(end-K+1:end));
end
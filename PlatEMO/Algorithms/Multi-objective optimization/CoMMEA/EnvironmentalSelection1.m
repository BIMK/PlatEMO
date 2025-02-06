function [Population,Fitness] = EnvironmentalSelection1(Population,N,Operation,EvoState)
% The environmental selection of CoMMEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% SPEA2 selection
    Fitness = CalFitness(Population,Operation);
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sortrows([Fitness' -Crowding(Population.decs)]);
        Population = Population(Rank(1:N));
    else
        Population = Population(Next);
        while length(Population)>N
            dist = sort(pdist2(Population.decs,Population.decs));
            CrowdDis = sum(dist(1:3,:));
            [~,index] = min(CrowdDis);
            Population(index) = [];
        end
    end
    if EvoState<0.5
        Fitness=CalFitness(Population,Operation)-Crowding(Population.decs)';
    else
        Fitness=-Crowding(Population.decs)';
    end
end
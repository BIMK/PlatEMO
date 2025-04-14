function [Population,Fitness] = EnvironmentalSelection(Population,N)
% The environmental selection for single- and multi-objective optimization

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if isscalar(Population(1).obj)  % Single-objective optimization
        Fitness    = FitnessSingle(Population);
        [~,rank]   = sort(Fitness);
        Population = Population(rank(1:N));
        Fitness    = Fitness(rank(1:N));
    else                            % Multi-objective optimization
        % Non-dominated sorting
        [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
        Next = FrontNo < MaxFNo;
        % Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(Population.objs,FrontNo);
        % Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:N-sum(Next)))) = true;
        % Population for next generation
        Population = Population(Next);
        Fitness    = FrontNo(Next);
    end
end
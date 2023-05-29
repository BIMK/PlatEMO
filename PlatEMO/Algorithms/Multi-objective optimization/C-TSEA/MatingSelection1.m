function MatingPool = MatingSelection1(Population,RefPoint,Range)
% The mating selection of AR-MOEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the fitness of each feasible solution based on IGD-NS
    % Calculate the distance between each solution and point
    N = length(Population);
    Distance    = CalDistance(Population.objs-repmat(Range(1,:),N,1),RefPoint);
    Convergence = min(Distance,[],2);
    [dis,rank]  = sort(Distance,1);
    % Calculate the fitness of noncontributing solutions
    Noncontributing = true(1,N);
    Noncontributing(rank(1,:)) = false;
    METRIC   = sum(dis(1,:)) + sum(Convergence(Noncontributing));
    fitness  = inf(1,N);
    fitness(Noncontributing) = METRIC - Convergence(Noncontributing);
    % Calculate the fitness of contributing solutions
    for p = find(~Noncontributing)
        temp = rank(1,:) == p;
        noncontributing = false(1,N);
        noncontributing(rank(2,temp)) = true;
        noncontributing = noncontributing & Noncontributing;
        fitness(p) = METRIC - sum(dis(1,temp)) + sum(dis(2,temp)) - sum(Convergence(noncontributing));
    end

    %% Combine the fitness of feasible solutions with the fitness of infeasible solutions
    Fitness = -inf(1,length(Population));
    
    %% Binary tournament selection
    MatingPool = TournamentSelection(2,length(Population),-Fitness);
end
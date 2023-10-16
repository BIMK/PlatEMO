function MatingPool = MatingSelection(Obj,RefPoint,Range)
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
        N = size(Obj,1);
        Distance    = CalDistance(Obj-repmat(Range(1,:),N,1),RefPoint);%N*NR
        Convergence = min(Distance,[],2);% The minimum Angle between each solution and the reference point
        [dis,rank]  = sort(Distance,1);% In columns, for reference points
        % Calculate the fitness of noncontributing solutions
        Noncontributing = true(1,N);
        Noncontributing(rank(1,:)) = false;
        METRIC   = sum(dis(1,:)) + sum(Convergence(Noncontributing));
        % The sum of the minimum angles of all reference points plus the sum of the minimum angles of all non-contributing solutions
        fitness  = inf(1,N);
        fitness(Noncontributing) = METRIC - Convergence(Noncontributing);
        % Calculate the fitness of contributing solutions
        for p = find(~Noncontributing)
            temp = rank(1,:) == p;
            noncontributing = false(1,N);
            noncontributing(rank(2,temp)) = true;% After removing this contribution point, the new contribution point Index
            noncontributing = noncontributing & Noncontributing;% Judge whether the new contribution point is the original non-contribution point or not
            fitness(p) = METRIC - sum(dis(1,temp)) + sum(dis(2,temp)) - sum(Convergence(noncontributing));
            % When this contribution point is removed, the METRIC adjusts
        end

    %% Combine the fitness of feasible solutions with the fitness of infeasible solutions
    Fitness = fitness;
    
    %% Binary tournament selection
    CV=zeros(N,1);
    MatingPool = TournamentSelection(2,ceil(N/2)*2,CV,-Fitness);
end
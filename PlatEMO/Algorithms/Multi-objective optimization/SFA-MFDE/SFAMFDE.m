classdef SFAMFDE < ALGORITHM
% <2025> <multi/many> <real/integer> <expensive>
% Scalarization function approximation based multifactorial differential evolution algorithm
% F     --- 0.5 --- The scaling factor for inner DE
% CR    --- 0.9 --- The crossover rate for inner DE
% omega ---  20 --- The maximum number of generations for MFDE-based search
% RMP   --- 0.3 --- Random Mating Probability
% nt    ---   2 --- The maximum number of aggregated subproblems
% nr    ---   2 --- The maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% Y. Horaguchi and M. Nakata. High-dimensional expensive multiobjective
% optimization using a surrogate-assisted multifactorial evolutionary
% algorithm. Proceedings of the Genetic and Evolutionary Computation
% Conference, 2025, 572-580.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yuma Horaguchi (email: horaguchi-yuma-hv@ynu.jp)

    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            [F, CR, omega, RMP, nt, nr] = Algorithm.ParameterSet(0.5, 0.9, 20, 0.3, 2, 2);

            %% Generate the weight vectors
            [W, Problem.N] = UniformPoint(Problem.N, Problem.M);

            %% Detect the neighbours of each solution
            T      = ceil(Problem.N / 10);
            B      = pdist2(W, W);
            [~, B] = sort(B, 2);
            B      = B(:, 1 : T);

            %% Initialize population
            PopDec     = UniformPoint(Problem.N, Problem.D, 'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper - Problem.lower, Problem.N, 1) .* PopDec + repmat(Problem.lower, Problem.N, 1));
            Arc        = Population;
            Z          = min(Population.objs, [], 1);

            %% Optimization
            while Algorithm.NotTerminated(Arc)
                for i = 1 : Problem.N
                    [~, UniqueID] = unique(Arc.decs, 'stable', 'rows');
                    ArcUnique     = Arc(UniqueID);

                    % Aggregate subproblem randomly
                    Remain = setdiff(1 : Problem.N, i);
                    TaskID = [i, Remain(randperm(Problem.N - 1, nt - 1))];

                    % Model Construction
                    RBFModel = {};
                    for j = 1 : length(TaskID)
                        ArcTch        = max(abs(ArcUnique.objs - repmat(Z, length(ArcUnique), 1)) .* W(TaskID(j), :), [], 2);
                        [tr_y, SrtID] = sort(ArcTch);
                        tr_x          = ArcUnique(SrtID(1 : Problem.N)).decs;
                        tr_y          = tr_y(1 : Problem.N);
                        pair          = pdist2(tr_x, tr_x);
                        D_max         = max(max(pair, [], 2));
                        spread        = D_max * (Problem.D * Problem.N) ^ (-1 / Problem.D);
                        RBFModel{j}   = newrbe(transpose(tr_x), transpose(tr_y), spread);
                    end

                    % MFDE-based search
                    P             = B(TaskID, randperm(end));
                    MTO           = MFDE(Problem, RBFModel, RMP, F, CR, nt);
                    Candidates    = MTO.run(Population(P).decs, omega);
                    [~, UniqueID] = unique(Candidates.decs, 'stable', 'rows');
                    Candidates    = Candidates(UniqueID);
                    [~, SrtID]    = sort(Candidates.objs, 1);
                    [~, BestID]   = min(sum(SrtID, 2));
                    OffDec        = Candidates(BestID).decs;

                    % Evaluate offspting
                    Offspring = Problem.Evaluation(OffDec);

                    % Update the reference point
                    Z = min(Z, Offspring.obj);

                    % Update population and archive
                    P     = unique(P, 'stable');
                    P     = P(randperm(end));
                    g_old = max(abs(Population(P).objs - repmat(Z, length(P), 1)) .* W(P, :), [], 2);
                    g_new = max(repmat(abs(Offspring.obj - Z), length(P), 1) .* W(P, :), [], 2);
                    Population(P(find(g_old >= g_new, nr))) = Offspring;
                    Arc   = [Arc, Offspring];

                    % Check termination criteria
                    Algorithm.NotTerminated(Arc);
                end
            end
        end
    end
end
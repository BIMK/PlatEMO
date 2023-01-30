classdef CLIA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Evolutionary algorithm with cascade clustering and reference point incremental learning

%------------------------------- Reference --------------------------------
% H. Ge, M. Zhao, L. Sun, Z. Wang, G. Tan, Q. Zhang, and C. L. P. Chen, A
% many-objective evolutionary algorithm with two interacting processes:
% Cascade clustering and reference point incremental learning, IEEE
% Transactions on Evolutionary Computation, 2019, 23(4): 572-586.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            global stable_threshold delta crowding_pick_flag;
            [stable_threshold, delta] = Algorithm.ParameterSet([0, 0, 0], 2 * Problem.N);
            [Z, P, A, S, SVM] = initialize(Problem);
            while Algorithm.NotTerminated(S)
                MatingPool = TournamentSelection(2, Problem.N, sum(max(0, P.cons), 2)); 
                Offspring  = OperatorGA(Problem,P(MatingPool));
                A = update_archive(A, [P, Offspring], Z, ceil(0.33 * Problem.M * Problem.N), Problem);
                % SELECTION OF INDIVIDUALS
                [P, ICA, ICN] = cascade_cluster([P, Offspring], Z, 'PDM', Problem.N, Problem.FE < Problem.maxFE);
                % ADAPTATION OF REFERENCE VECTORS     
                [Z, SVM] = incremental_learn(Z, ICA, ICN, A, SVM, Problem);
                if Problem.FE >= Problem.maxFE && crowding_pick_flag
                    S = crowding_pick(update_archive(A, P, Z, [], Problem), Problem.N, 'precise');
                else
                    S = P;
                end
            end
        end
    end
end
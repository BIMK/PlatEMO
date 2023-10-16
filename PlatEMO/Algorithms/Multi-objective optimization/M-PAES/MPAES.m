classdef MPAES < ALGORITHM
% <multi> <real/integer>
% Memetic algorithm with Pareto archived evolution strategy
% l_fails   ---  5 --- Maximum number of consecutive failing local moves
% l_opt     --- 10 --- Maximum number of local moves
% cr_trials --- 20 --- Number of crossover trials
% div       --- 10 --- The number of divisions in each objective

%------------------------------- Reference --------------------------------
% J. D. Knowles and D. W. Corne, M-PAES: A memetic algorithm for
% multiobjective optimization, Proceedings of the IEEE Congress on
% Evolutionary Computation, 2000, 325-332.
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
            %% Parameter setting
            [l_fails,l_opt,cr_trials,div] = Algorithm.ParameterSet(5,10,20,10);

            %% Generate random population
            P = Problem.Initialization();
            G = P(NDSort(P.objs,1)==1);

            %% Optimization
            while Algorithm.NotTerminated(G)
                for i = 1 : Problem.N
                    H = G(~all(G.objs<=repmat(P(i).obj,length(G),1),2));
                    [P(i),G] = PAES(Problem,P(i),G,[H,P(i)],Problem.N,l_fails,l_opt,div);
                end
                P1(1:Problem.N) = SOLUTION();
                for i = 1 : Problem.N
                    for r = 1 : cr_trials
                        Combine = [P,G];
                        parents = Combine(randperm(length(Combine),2));
                        c       = OperatorGAhalf(Problem,parents,{1,20,0,0});
                        [G,dominated,GCrowd,cCrowd,pCrowd] = UpdateArchive(G,c,parents,Problem.N,div);
                        if ~dominated && any(cCrowd<=pCrowd)
                            break;
                        end
                    end
                    if dominated
                        c = G(TournamentSelection(2,1,GCrowd));
                    end
                    P1(i) = c;
                end
                P = P1;
            end
        end
    end
end
classdef PRCEA < ALGORITHM
% <2026> <multi/many> <real/integer> <large/none> <constrained/none>
% Promising region-guided large-scale constrained MOEA

%------------------------------- Reference --------------------------------
% X. Zhong, X. Yao, K. Qiao, D. Gong, and Y. Jin. A large-scale constrained
% multi-objective evolutionary algorithm with promising region detection
% and diversity enhancement. IEEE Transactions on Evolutionary Computation,
% 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Pop{1} = Problem.Initialization();
            Pop{2} = Problem.Initialization();
            Off    = {[],[]};
            stage            = 1;
            lgap             = 20;
            change_threshold = 1e-3;
            NP               = Problem.N;
            cnt              = 0;
    
            %% Optimization
            while Algorithm.NotTerminated(Pop{1})
                cnt = cnt + 1;
                % Population Update
                [Pop{1},Fit{1}] = CDP_EPD_Update([Pop{1}, Off{1},Off{2}],NP, true);
                if stage == 1
                    [Pop{2},Fit{2}] = CDP_EPD_Update([Pop{2}, Off{1},Off{2}],NP,false);
                else
                    [Pop{2},Fit{2}] = PRDD_Update(Pop{1},[Pop{2}, Off{1},Off{2}],NP);
                end
                % Udate ideal point
                for i = 1 : 2
                    zmin{i} = min([Pop{i}.objs], [], 1) - 1e-6;
                end
                % Update stage
                if stage == 1
                    Objs = Pop{1}.objs;
                    ideal_points(cnt,:) = min(Objs, [], 1);
                    nadir_points(cnt,:) = max(Objs, [], 1);
                    mean_points(cnt,:)  = mean(Objs, 1);
                    if cnt >= lgap
                        max_change = calc_maxchange(ideal_points, nadir_points, mean_points, cnt, lgap);
                        if max_change <= change_threshold
                            stage = 2;
                        end
                    end
                end
                % Offspring generation
                Off = cell(1, 2);
                for i = 1 : 2
                    if stage == 1
                        Off{i} = AES(Problem, Pop{i}, Fit{i}, zmin{i});
                    else
                        MatingPool          = [Pop{i}(randsample(NP,NP))];
                        [Mate1,Mate2,Mate3] = Neighbor_Pairing_Strategy(MatingPool,zmin{i});
                        Off{i}(1:NP/2)      = OperatorDE_rand_1(Problem,Mate1(1:NP/2), Mate2(1:NP/2), Mate3(1:NP/2));
                        Off{i}(1+NP/2:NP)   = OperatorDE_pbest_1_main(Problem, Pop{i}, NP/2, Fit{i}, 0.1);
                    end
                end
                if Problem.FE >= Problem.maxFE
                    [Pop{1},~] = CDP_EPD_Update([Pop{1}, Off{1},Off{2}],Problem.N, true);
                end
            end
        end
    end
end
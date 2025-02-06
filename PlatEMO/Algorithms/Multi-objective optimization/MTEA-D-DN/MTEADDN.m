classdef MTEADDN < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <multitask>
% Multiobjective multitask evolutionary algorithm based on decomposition with dual neighborhoods

%------------------------------- Reference --------------------------------
% X. Wang, Z. Dong, L. Tang, and Q. Zhang. Multiobjective multitask
% optimization-neighborhood as a bridge for knowledge transfer. IEEE
% Transactions on Evolutionary Computation, 2023, 27(1): 155-169.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            Beta  = 0.2;
            F     = 0.5;
            CR    = 0.9;
            MuM   = 20;
            ProbT = size(Problem.SubM,2);
            % Initialize
            for t = 1 : ProbT
                % Generate the weight vectors
                [W{t}, N{t}] = UniformPoint(Problem.N/2, Problem.SubM(t));
                DT{t} = ceil(N{t} / 20);

                % Detect the neighbours of each solution
                B{t}      = pdist2(W{t}, W{t});
                [~, B{t}] = sort(B{t}, 2);
                B{t}      = B{t}(:, 1:DT{t});
                Dec       = [rand(Problem.N/2, max(Problem.SubD)),t*ones(Problem.N/2,1)];
                SubPopulation{t} = Problem.Evaluation(Dec);           
                Z{t} = min(SubPopulation{t}.objs, [], 1);
            end
            for t = 1 : ProbT
                % Detect the second neighbours of each solution
                tar_pool = 1:ProbT; tar_pool(t) = [];
                for i = 1 : size(B{t}, 1)
                    B2k{t}(i) = tar_pool(randi(ProbT - 1));
                    B2{t, i}  = randperm(size(W{B2k{t}(i)}, 1), DT{t});
                end
            end

            %% Optimization            
            while Algorithm.NotTerminated([SubPopulation{:}])
                % Generation
                for t = 1 : ProbT
                    for i = 1 : N{t}
                        % Choose the parents
                        if rand() < Beta
                            P1    = B{t}(i, :);
                            P2    = B2{t, i};
                            tasks = [t * ones(1, length(P1)), B2k{t}(i) * ones(1, length(P2))];
                            P     = [P1, P2];
                            rndpm = randperm(length(tasks));
                            tasks = tasks(rndpm);
                            P     = P(rndpm);

                            % Generate an offspring
                            offspringDec = OperatorDE(Problem,SubPopulation{t}(i).dec,SubPopulation{tasks(1)}(P(1)).dec,SubPopulation{tasks(2)}(P(2)).dec,{CR,F,1,MuM});
                            offspringDec(end) = t;
    
                            if rand() < 0.5
                                % Knowledge transfer
                                k = B2k{t}(i);
                                offspringDec(end) = k;
                                offspring = Problem.Evaluation(offspringDec);
                                Z{k}      = min(Z{k}, offspring.obj);
                                g_old     = max(abs(SubPopulation{k}(P2).objs - repmat(Z{k}, length(P2), 1)) .* W{k}(P2, :), [], 2);
                                g_new     = max(repmat(abs(offspring.obj - Z{k}), length(P2), 1) .* W{k}(P2, :), [], 2);
                                SubPopulation{k}(P2(g_old >= g_new)) = offspring;
    
                                if all(g_old < g_new)
                                    tar_pool  = 1:ProbT; tar_pool(t) = [];
                                    B2k{t}(i) = tar_pool(randi(ProbT - 1));
                                    B2{t, i}  = randperm(size(W{B2k{t}(i)}, 1), DT{t});
                                elseif any (g_old >= g_new)
                                    B2{t, i} = B2{t, i}(g_old >= g_new);
                                end
                            else
                                offspring = Problem.Evaluation(offspringDec);
                                Z{t}      = min(Z{t}, offspring.obj);
                                g_old     = max(abs(SubPopulation{t}(P1).objs - repmat(Z{t}, length(P1), 1)) .* W{t}(P1, :), [], 2);
                                g_new     = max(repmat(abs(offspring.obj - Z{t}), length(P1), 1) .* W{t}(P1, :), [], 2);
                                SubPopulation{t}(P1(g_old >= g_new)) = offspring;
                            end
                        else
                            P = randperm(N{t});
                            offspringDec = OperatorDE(Problem,SubPopulation{t}(i).dec,SubPopulation{t}(P(1)).dec,SubPopulation{t}(P(2)).dec,{CR,F,1,MuM});
                            offspringDec(end) = t;
                            offspring    = Problem.Evaluation(offspringDec);
                            Z{t}         = min(Z{t}, offspring.obj);
                            g_old        = max(abs(SubPopulation{t}(P).objs - repmat(Z{t}, length(P), 1)) .* W{t}(P, :), [], 2);
                            g_new        = max(repmat(abs(offspring.obj - Z{t}), length(P), 1) .* W{t}(P, :), [], 2);
                            SubPopulation{t}(P(g_old >= g_new)) = offspring;
                        end
                    end
                end
            end
        end
    end
end
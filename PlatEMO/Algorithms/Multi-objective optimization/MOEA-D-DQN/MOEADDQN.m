classdef MOEADDQN < ALGORITHM
% <multi/many> <real/integer>
% MOEA/D based on deep Q-network

%------------------------------- Reference --------------------------------
% Y. Tian, X. Li, H. Ma, X. Zhang, K. C. Tan, and Y. Jin, Deep
% reinforcement learning based adaptive operator selection for evolutionary
% multi-objective optimization, IEEE Transactions on Emerging Topics in
% Computational Intelligence, 2023, 7(4): 1051-1064.
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
            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Operators
            crossOp = RecRL(Problem, W, Problem.maxFE, Problem.N);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);
            % Utility for each subproblem
            Pi     = ones(Problem.N,1);
            oldObj = max(abs((Population.objs-Z)./W),[],2);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                for i = 1 : 5
                    % Search for candidates
                    Boundary  = find(sum(W<1e-3,2)==Problem.M-1)';
                    Candidate = [Boundary,TournamentSelection(10,floor(Problem.N/5)-length(Boundary),-Pi)];
                    Candidate = unique(Candidate);
                    
                    % For each candidate
                    for j = 1 : size(Candidate,2)
                        % Choose neighbour, or treat all as neighbour
                        % indices of neighbor is stored in P
                        if rand < 0.9
                            P = B(Candidate(j),randperm(end));
                        else
                            P = randperm(Problem.N);
                        end
                        % Choose offspring
                        [crossOp, Offspring] = crossOp.do(Population.decs, Candidate(j), P, Problem.FE);
                        % Mutate
                        Offspring = PolyMut(Offspring, Problem.lower, Problem.upper);
                        % Evaluate
                        Offspring = Problem.Evaluation(Offspring);
                        % Determine which to update
                        Z = min(Z,Offspring.obj);
                        g_old = max(abs((Population(P).objs-Z)./W(P,:)),[],2);
                        g_new = max(abs((Offspring.obj-Z)./W(P,:)),[],2);
                        update_idxs = (g_old>=g_new);
                        % Update and train cross operator
                        if sum(update_idxs) >= 1
                            FIR = (g_old(update_idxs) - g_new(update_idxs)) ./ g_old(update_idxs);
                            Population(P(g_old>=g_new)) = Offspring;
                            crossOp = crossOp.learn(sum(FIR));
                        end
                    end
                end
                if ~mod(ceil(Problem.FE/Problem.N),10)
                    % Update Pi for each solution
                    newObj    = max(abs((Population.objs-Z)./W),[],2);
                    DELTA     = (oldObj-newObj)./oldObj;
                    Temp      = DELTA < 0.001;
                    Pi(~Temp) = 1;
                    Pi(Temp)  = (0.95+0.05*DELTA(Temp)/0.001).*Pi(Temp);
                    oldObj    = newObj;
                end
            end
        end
    end
end
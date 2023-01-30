classdef ENSMOEAD < ALGORITHM
% <multi/many> <real/integer>
% Ensemble of different neighborhood sizes based MOEA/D
% NS --- 25:25:100 --- Set of neighborhood sizes
% LP ---        50 --- Learning period

%------------------------------- Reference --------------------------------
% S. Zhao, P. N. Suganthan, and Q. Zhang, Decomposition-based multi-
% objective evolutionary algorithm with an ensemble of neighborhood sizes,
% IEEE Transactions on Evolutionary Computation, 2012, 16(3): 442-446.
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
            [NS,LP] = Algorithm.ParameterSet(25:25:100,50);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            % Maximum number of solutions replaced by each offspring
            nr = ceil(Problem.N/100);

            %% Detect all the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);

            %% Generate random population
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);
            % Utility for each subproblem
            Pi = ones(Problem.N,1);
            % Old Tchebycheff function value of each solution on its subproblem
            oldObj = max(abs((Population.objs-repmat(Z,Problem.N,1)).*W),[],2);

            %% Optimization
            p           = ones(1,length(NS))./length(NS);
            FEs         = zeros(1,length(NS));
            FEs_success = zeros(1,length(NS));
            while Algorithm.NotTerminated(Population)
                % Select neighborhood size for each subproblem
                ns = RouletteWheelSelection(Problem.N,1./p);

                % Apply MOEA/D-DRA for one generation
                for subgeneration = 1 : 5
                    % Choose I
                    Bounday = find(sum(W<1e-3,2)==Problem.M-1)';
                    I = [Bounday,TournamentSelection(10,floor(Problem.N/5)-length(Bounday),-Pi)];

                    % For each solution in I
                    for i = I
                        % Choose the parents
                        if rand < 0.9
                            P = B(i,randperm(NS(ns(i))));
                        else
                            P = randperm(Problem.N);
                        end

                        % Generate an offspring
                        Offspring = OperatorDE(Problem,Population(i),Population(P(1)),Population(P(2)));

                        % Update the ideal point
                        Z = min(Z,Offspring.obj);

                        % Update the solutions in P by Tchebycheff approach
                        g_old   = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                        g_new   = max(repmat(abs(Offspring.obj-Z),length(P),1).*W(P,:),[],2);
                        replace = find(g_old>=g_new,nr);
                        Population(P(replace)) = Offspring;
                        if ~isempty(replace)
                            FEs_success(ns(i)) = FEs_success(ns(i)) + 1;
                        end
                        FEs(ns(i)) = FEs(ns(i)) + 1;
                    end
                end
                if ~mod(ceil(Problem.FE/Problem.N),10)
                    % Update Pi for each solution
                    newObj    = max(abs((Population.objs-repmat(Z,Problem.N,1)).*W),[],2);
                    DELTA     = (oldObj-newObj)./oldObj;
                    Temp      = DELTA < 0.001;
                    Pi(~Temp) = 1;
                    Pi(Temp)  = (0.95+0.05*DELTA(Temp)/0.001).*Pi(Temp);
                    oldObj    = newObj;
                end
                if ~mod(ceil(Problem.FE/Problem.N),LP)
                    % Update the probability of choosing each neighborhood size
                    R           = FEs_success./FEs;
                    p           = R./sum(R);
                    FEs         = zeros(1,length(NS));
                    FEs_success = zeros(1,length(NS));
                end
            end
        end
    end
end
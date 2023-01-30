classdef MOEADFRRMAB < ALGORITHM
% <multi/many> <real/integer>
% MOEA/D with fitness-rate-rank-based multiarmed bandit
% C --- 5 --- Scaling factor in bandit-based operator selection
% W ---   --- Size of sliding window
% D --- 1 --- Decaying factor in calculating credit value

%------------------------------- Reference --------------------------------
% K. Li, A. Fialho, S. Kwong, and Q. Zhang, Adaptive operator selection
% with bandits for a multiobjective evolutionary algorithm based on
% decomposition, IEEE Transactions on Evolutionary Computation, 2014,
% 18(1): 114-130.
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
            [C,W,D] = Algorithm.ParameterSet(5,ceil(Problem.N/2),1);

            %% Generate the weight vectors
            [Weight,Problem.N] = UniformPoint(Problem.N,Problem.M);
            % Size of neighborhood
            T = 20;
            % Maximum number of solutions replaced by each offspring
            nr = 2;

            %% Detect the neighbours of each solution
            B = pdist2(Weight,Weight);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);
            % Utility for each subproblem
            Pi = ones(Problem.N,1);
            % Old Tchebycheff function value of each solution on its subproblem
            oldObj = max(abs((Population.objs-repmat(Z,Problem.N,1)).*Weight),[],2);

            %% Optimization
            FRR = zeros(1,4);	% Credit value of each operator
            SW  = zeros(2,W);	% Sliding window
            while Algorithm.NotTerminated(Population)
                for subgeneration = 1 : 5
                    % Choose I
                    Bounday = find(sum(Weight<1e-3,2)==Problem.M-1)';
                    I = [Bounday,TournamentSelection(10,floor(Problem.N/5)-length(Bounday),-Pi)];

                    % For each solution in I
                    for i = I
                        % Bandit-based operator selection
                        op = FRRMAB(FRR,SW,C);

                        % Choose the parents
                        if rand < 0.9
                            P = B(i,randperm(end));
                        else
                            P = randperm(Problem.N);
                        end

                        % Generate an offspring
                        Offspring = FourDE(Problem,op,Population(i),Population(P(1)),Population(P(2)),Population(P(3)),Population(P(4)),Population(P(5)));

                        % Update the ideal point
                        Z = min(Z,Offspring.obj);

                        % Update the solutions in P by Tchebycheff approach
                        g_old   = max(abs(Population(P).objs-repmat(Z,length(P),1)).*Weight(P,:),[],2);
                        g_new   = max(repmat(abs(Offspring.obj-Z),length(P),1).*Weight(P,:),[],2);
                        replace = find(g_old>=g_new,nr);
                        Population(P(replace)) = Offspring;
                        FIR = sum((g_old(replace)-g_new(replace))./g_old(replace));
                        SW  = [SW(1,2:end),op;SW(2,2:end),FIR];
                        FRR = CreditAssignment(SW,D);
                    end
                end
                if ~mod(ceil(Problem.FE/Problem.N),10)
                    % Update Pi for each solution
                    newObj    = max(abs((Population.objs-repmat(Z,Problem.N,1)).*Weight),[],2);
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
classdef MOEADDRA < ALGORITHM
% <multi/many> <real/integer>
% MOEA/D with dynamical resource allocation

%------------------------------- Reference --------------------------------
% Q. Zhang, W. Liu, and H. Li, The performance of a new version of MOEA/D
% on CEC09 unconstrained MOP test instances, Proceedings of the IEEE
% Congress on Evolutionary Computation, 2009, 203-208.
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
            % Size of neighborhood
            T  = ceil(Problem.N/10);
            % Maximum number of solutions replaced by each offspring
            nr = ceil(Problem.N/100);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);
            % Utility for each subproblem
            Pi = ones(Problem.N,1);
            % Old Tchebycheff function value of each solution on its subproblem
            oldObj = max(abs((Population.objs-repmat(Z,Problem.N,1)).*W),[],2);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                for subgeneration = 1 : 5
                    % Choose I
                    Bounday = find(sum(W<1e-3,2)==Problem.M-1)';
                    I = [Bounday,TournamentSelection(10,floor(Problem.N/5)-length(Bounday),-Pi)];

                    % For each solution in I
                    for i = I
                        % Choose the parents
                        if rand < 0.9
                            P = B(i,randperm(end));
                        else
                            P = randperm(Problem.N);
                        end

                        % Generate an offspring
                        Offspring = OperatorDE(Problem,Population(i),Population(P(1)),Population(P(2)));

                        % Update the ideal point
                        Z = min(Z,Offspring.obj);

                        % Update the solutions in P by Tchebycheff approach
                        g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                        g_new = max(repmat(abs(Offspring.obj-Z),length(P),1).*W(P,:),[],2);
                        Population(P(find(g_old>=g_new,nr))) = Offspring;
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
            end
        end
    end
end
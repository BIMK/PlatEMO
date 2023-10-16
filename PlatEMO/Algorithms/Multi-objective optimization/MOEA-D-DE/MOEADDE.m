classdef MOEADDE < ALGORITHM
% <multi/many> <real/integer>
% MOEA/D based on differential evolution
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% H. Li and Q. Zhang, Multiobjective optimization problems with complicated
% Pareto sets, MOEA/D and NSGA-II, IEEE Transactions on Evolutionary
% Computation, 2009, 13(2): 284-302.
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
            [delta,nr] = Algorithm.ParameterSet(0.9,2);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % For each solution
                for i = 1 : Problem.N
                    % Choose the parents
                    if rand < delta
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
        end
    end
end
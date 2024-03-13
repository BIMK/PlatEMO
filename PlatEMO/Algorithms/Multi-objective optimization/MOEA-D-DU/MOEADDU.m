classdef MOEADDU < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% MOEA/D with a distance based updating strategy
% delta --- 0.9 --- The probability of choosing parents locally
% K     ---   5 --- Number of nearest weight vectors

%------------------------------- Reference --------------------------------
% Y. Yuan, H. Xu, B. Wang, B. Zhang, and X. Yao, Balancing convergence and
% diversity in decomposition-based many-objective optimizers, IEEE
% Transactions on Evolutionary Computation, 2016, 20(2): 180-198.
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
            [delta,K] = Algorithm.ParameterSet(0.9,5);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            % Ideal and nadir points
            z    = min(Population.objs,[],1);
            znad = max(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Normalization
                [~,z,znad] = Normalization(Population.objs,z,znad);
                % For each solution
                for i = 1 : Problem.N
                    % Choose the parent
                    if rand < delta
                        P = B(i,randi(size(B,2)));
                    else
                        P = randi(Problem.N);
                    end

                    % Generate an offspring
                    Offspring = OperatorGAhalf(Problem,Population([i,P]));

                    % Find the K nearest weight vectors of the offspring
                    [~,rank] = sort(1-pdist2(Offspring.obj,W,'cosine'),'descend');
                    P        = rank(1:K);

                    % Update the K nearest parents by modified Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(z,length(P),1))./repmat(znad-z,K,1)./W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-z)./(znad-z),length(P),1)./W(P,:),[],2);
                    Population(P(find(g_old>=g_new,1))) = Offspring;
                end
            end
        end
    end
end
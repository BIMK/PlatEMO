classdef MOEADVOV < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% MOEA/D with virtual objective vectors
% G     ---  100 --- The weight update generation
% C     ---    9 --- The weight update count
% theta --- 0.02 --- Threshold

%------------------------------- Reference --------------------------------
% T. Takagi, K. Takadama, and H. Sato, Weight vector arrangement using
% virtual objective vectors in decomposition-based MOEA, Proceedings of
% the IEEE Congress on Evolutionary Computation, 2021, 1462-1469.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [G,C,theta] = Algorithm.ParameterSet(100,9,0.02);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M,'ILD');
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);

             %% Set bifurcation point and external population
            BP = Problem.N * G * C;
            EP = Population;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % For each solution
                Offsprings(1:Problem.N) = SOLUTION();
                for i = 1 : Problem.N
                    % Choose the parents
                    P = B(i,randperm(size(B,2)));

                    % Generate an offspring
                    Offsprings(i) = OperatorGAhalf(Problem,Population(P(1:2)));

                    % Update the ideal point
                    Z = min(Z,Offsprings(i).obj);

                    % Update the solutions in P by Modified Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                    g_new = max(repmat(abs(Offsprings(i).obj-Z),T,1)./W(P,:),[],2);
                    Population(P(g_old>=g_new)) = Offsprings(i);
                end

                if Problem.FE <= BP
                    EP = [EP,Offsprings];
                    EP = EP(NDSort(EP.objs,1)==1);
                    if length(EP) > 5000
                        EP(1:length(EP)-5000) = [];
                    end
                    if ~mod(ceil(Problem.FE/Problem.N),G)
                        % Generate the virtual objective vectors
                        VOV = generateVOV(EP.objs,theta);
                        % Update the weight vectors
                        [W,B] = updateWeight(EP.objs,VOV,Problem.N);
                        % Update the population 
                        obj = abs(EP.objs-repmat(Z,length(EP),1));
                        for i = 1 : Problem.N
                            [~,I] = min(max(obj./W(i,:),[],2));
                            Population(i) = EP(I);
                        end
                    end
                end
            end
        end
    end
end 
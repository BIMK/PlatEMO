classdef MOEADDCWV < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% MOEA/D with distribution control of weight vector set
% p --- -1 --- The intermediate objective value

%------------------------------- Reference --------------------------------
% T. Takagi, K. Takadama, and H. Sato, A distribution control of weight 
% vector set for multi-objective evolutionary algorithms, Proceedings of 
% the Bio-inspired Information and Communication Technologies, Lecture 
% Notes of the Institute for Computer Sciences, Social Informatics and 
% Telecommunications Engineering, 2019, 70-80.
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
            p = Algorithm.ParameterSet(-1);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Set the weight vectors
            [W,W0] = setWeight(W,p);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if ~isempty(W0)
                    W = updateWeight(Population.objs, W0);
                end
                % For each solution
                for i = 1 : Problem.N
                    % Choose the parents
                    P = B(i,randperm(size(B,2)));

                    % Generate an offspring
                    Offspring = OperatorGAhalf(Problem,Population(P(1:2)));

                    % Update the ideal and nadir point
                    Z = min(Z,Offspring.obj);
                    Zmax  = max(Population.objs,[],1);

                    % Update the solutions in P by Modified Tchebycheff approach with normalization
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1)./W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1)./W(P,:),[],2);
                    Population(P(g_old>=g_new)) = Offspring;
                end
            end
        end
    end
end 
classdef WVMOEAP < ALGORITHM
% <multi> <real/integer>
% Weight vector based multi-objective optimization algorithm with preference
% Points ---      --- Set of preferred points
% b      --- 0.05 --- Range of preferred region

%------------------------------- Reference --------------------------------
% X. Zhang, X. Jiang, and L. Zhang, A weight vector based multi-objective
% optimization algorithm with preference, Acta Electronica Sinica
% (Chinese), 2016, 44(11): 2639-2645.
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
            [Points,b] = Algorithm.ParameterSet(ones(1,Problem.M),0.05);

            %% Generate the weight vectors and random population
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Map the weight vectors
            Dis   = pdist2(W,W);
            B     = zeros(Problem.N,T);
            Group = ceil((1:Problem.N)/Problem.N*size(Points,1));
            for i = unique(Group)
                % The weight vectors around the i-th preferred point
                Current = find(Group==i);
                % Map the weight vectors
                W(Current,:) = 2*b.*W(Current,:) + Points(i,:) + b;
                % Detect the neighbours of each vector
                [~,rank]     = sort(Dis(Current,Current),2);
                B(Current,:) = Current(rank(:,1:T));
            end

            %% Generate random population
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % For each group
                for i = unique(Group)
                    Current = find(Group==i);
                    % Generate an offspring for each solution
                    P = zeros(length(Current),2);
                    for j = 1 : size(P,1)
                        if rand < 0.9
                            P(j,:) = B(j,randperm(size(B,2),2));
                        else
                            P(j,:) = Current(randperm(length(Current),2));
                        end
                    end
                    Offspring = OperatorDE(Problem,Population(Current),Population(P(:,1)),Population(P(:,2)));
                    % Environmental selection
                    Z = min(Z,min(Offspring.objs,[],1));
                    Population(Current) = STM([Population(Current),Offspring],W(Current,:),Z,Z+ones(1,Problem.M));
                end
            end
        end
    end
end
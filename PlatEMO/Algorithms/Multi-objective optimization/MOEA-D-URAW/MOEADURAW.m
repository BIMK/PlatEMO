classdef MOEADURAW < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% MOEA/D with uniform randomly adaptive weights
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% L. R. C. Farias and A. F. R. Araujo, Many-objective evolutionary
% algorithm based on decomposition with random and adaptive weights. In
% Proceedings of the IEEE International Conference on Systems, Mans and
% Cybernetics, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [delta,nr] = Algorithm.ParameterSet(0.9,2);

            %% Generate the weight vectors
            [W,Problem.N] = UniformlyRandomlyPoint(Problem.N,Problem.M);
            % Transformation on W
            W = 1./W./repmat(sum(1./W,2),1,size(W,2));
            % Size of neighborhood
            T = ceil(Problem.N/10);
            % Size of external elite
            nEP = ceil(Problem.N*2);
            % Ratio of updated weight vectors
            nus = 0.05;

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);

            %% Optimization
            EP = [];
            adaptation_moment = round(ceil(Problem.maxFE/Problem.N)*0.05);
            while Algorithm.NotTerminated(Population)
                % For each solution	
                Offsprings(1:Problem.N) = SOLUTION();
                for i = 1 : Problem.N
                    % Choose the parents
                    if rand < delta
                        P = B(i,randperm(size(B,2)));
                    else
                        P = randperm(Problem.N);
                    end

                    % Generate an offspring
                    Offsprings(i) = OperatorGAhalf(Problem,Population(P(1:2)));

                    % Update the ideal point
                    Z = min(Z,Offsprings(i).obj);

                    % Update the solutions in P by Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offsprings(i).obj-Z),length(P),1).*W(P,:),[],2);
                    Population(P(find(g_old>=g_new,nr))) = Offsprings(i);
                end
                if Problem.FE/Problem.maxFE <= 0.9
                    if isempty(EP)
                        EP = updateEP(Population,Offsprings,nEP);
                    else
                        EP = updateEP(EP,Offsprings,nEP);
                    end
                    if mod(ceil(Problem.FE/Problem.N),adaptation_moment) == 0
                        % Adaptive weight adjustment          
                        [Population,W] = updateWeight(Population,W,Z,EP,nus*Problem.N); 
                        B = pdist2(W,W);
                        [~,B] = sort(B,2);
                        B = B(:,1:T);
                    end
                end
            end
        end
    end
end
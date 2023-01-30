classdef SMEA < ALGORITHM
% <multi> <real/integer>
% Self-organizing multiobjective evolutionary algorithm
% D    ---     --- Number of neurons in each dimension of the latent space
% tau0 --- 0.7 --- Initial learning rate
% H    ---   5 --- Size of neighborhood mating pools

%------------------------------- Reference --------------------------------
% H. Zhang, A. Zhou, S. Song, Q. Zhang, X. Gao, and J. Zhang, A self-
% organizing multiobjective evolutionary algorithm, IEEE Transactions on
% Evolutionary Computation, 2016, 20(5): 792-806.
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
            [D,tau0,H] = Algorithm.ParameterSet(repmat(ceil(Problem.N.^(1/(Problem.M-1))),1,Problem.M-1),0.7,5);
            Problem.N  = prod(D);
            sigma0     = sqrt(sum(D.^2)/(Problem.M-1))/2;

            %% Generate random population
            Population = Problem.Initialization();
            FrontNo    = NDSort(Population.objs,inf);

            %% Initialize the SOM
            % Training set
            S = Population.decs;
            % Weight vector of each neuron
            W = S;
            % Position of each neuron
            D = arrayfun(@(S)1:S,D,'UniformOutput',false);
            eval(sprintf('[%s]=ndgrid(D{:});',sprintf('c%d,',1:length(D))))
            eval(sprintf('Z=[%s];',sprintf('c%d(:),',1:length(D))))
            % Distance between each two neurons in latent space
            LDis = pdist2(Z,Z);
            % H nearest neurons of each neuron in latent space
            [~,B] = sort(LDis,2);
            B     = B(:,2:min(H+1,end));

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Update SOM
                for s = 1 : size(S,1)
                    sigma  = sigma0*(1-(Problem.FE+s)/Problem.maxFE);
                    tau    = tau0*(1-(Problem.FE+s)/Problem.maxFE);
                    [~,u1] = min(pdist2(S(s,:),W));
                    U      = LDis(u1,:) < sigma;
                    W(U,:) = W(U,:) + tau.*repmat(exp(-LDis(u1,U))',1,size(W,2)).*(repmat(S(s,:),sum(U),1)-W(U,:));
                end

                % Associate each solution with a neuron
                A  = 1 : Problem.N;
                U  = 1 : Problem.N;
                XU = zeros(1,Problem.N);
                for i = 1 : Problem.N
                    x        = randi(length(A));
                    [~,u]    = min(pdist2(Population(A(x)).dec,W(U,:)));
                    XU(U(u)) = A(x);
                    A(x)     = [];
                    U(u)     = [];
                end

                % Evolution
                A = Population.decs;
                for u = 1 : Problem.N
                    drawnow('limitrate');
                    if rand < 0.9
                        Q = XU(B(u,:));
                    else
                        Q = 1 : Problem.N;
                    end
                    Q = Q(randperm(end,2));
                    y = OperatorDE(Problem,Population(u),Population(Q(1)),Population(Q(2)));
                    [Population,FrontNo] = Select(Population,FrontNo,y);
                end

                % Update the training set
                S = setdiff(Population.decs,A,'rows');
            end
        end
    end
end
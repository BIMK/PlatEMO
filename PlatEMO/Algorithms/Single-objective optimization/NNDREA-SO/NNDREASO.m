classdef NNDREASO < ALGORITHM
% <2024> <single> <binary> <large/none> <constrained/none> <sparse/none>
% Evolutionary algorithm with neural network-based dimensionality reduction
% lower ---  -1 --- Lower bound of network weights
% upper ---   1 --- Upper bound of network weights
% delta --- 0.5 --- Proportion of the first stage

%------------------------------- Reference --------------------------------
% Y. Tian, L. Wang, S. Yang, J. Ding, Y. Jin, and X. Zhang. Neural
% network-based dimensionality reduction for large-scale binary
% optimization with millions of variables. IEEE Transactions on
% Evolutionary Computation, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [lower, upper, delta] = Algorithm.ParameterSet(-1, 1, 0.5);
            % Extracting information from problems
            switch class(Problem)
                case 'KP'
                    Type = 1;
                    Instance = [Problem.P; Problem.W].';
                case 'MaxCut'
                    Type = 2;
                    AdjMat = Problem.Adj_mat;
                otherwise
                    error('The %s problem cannot be solved by this algorithm',class(Problem));
            end
            % Feature extraction
            if Type == 1       % Random perturbation for non-graph problems
                Instance = normrnd(Instance,0.1);
            elseif Type == 2   % Feature extraction for graph problems
                AdjMat   = sparse(AdjMat);
                D_mat    = zeros(size(AdjMat));
                D_mat(logical(eye(size(D_mat)))) = sum(AdjMat,2);
                D_mat    = D_mat ^ -0.5;
                D_mat    = sparse(D_mat);
                [V,DM]   = eig(full(D_mat*AdjMat*D_mat));
                [~,ind]  = sort(diag(DM));
                Instance = V(:,ind(1:10));
            end
            % Set neural network structures
            structure = [size(Instance, 2), 4, 1];
            s_list    = ones((size(structure, 2)-1)*2, 2) * -1;
            % Get list of neural network structures and dimension after
            % neural network weight flattening
            Dim = 0;
            for i = 1 : size(structure, 2)-1
                s_list(2*i-1,1) = structure(i);
                s_list(2*i-1,2) = structure(i+1);
                s_list(2*i,1)   = structure(i+1);
                Dim = Dim + structure(i) * structure(i+1) + structure(i+1);
            end
            % Get lower and upper bound of search space (all)
            lower = repmat(lower,1,Dim);
            upper = repmat(upper,1,Dim);

            %% Generate random population
            % Initialize population weights
            PopWeight = unifrnd(repmat(lower,Problem.N,1),repmat(upper,Problem.N,1));
            % Fully connected neural network for forward propagation
            Output = FCNForward(PopWeight, Instance, s_list);
            % Evaluate the population
            Population = Problem.Evaluation(Output);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FitnessSingle(Population));
                if Problem.FE <= Problem.maxFE * delta
                    % First stage
                    OffWeight  = OperatorReal(PopWeight(MatingPool, :),lower, upper);
                    Output     = FCNForward(OffWeight, Instance, s_list);
                    Offspring  = Problem.Evaluation(Output);
                    Population = [Population,Offspring];
                    PopWeight  = [PopWeight;OffWeight];
                    [~,rank]   = sort(FitnessSingle(Population));
                    Population = Population(rank(1:Problem.N));
                    PopWeight  = PopWeight(rank(1:Problem.N),:);
                else
                    % Second stage
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                    Population = [Population,Offspring];
                    [~,rank]   = sort(FitnessSingle(Population));
                    Population = Population(rank(1:Problem.N));
                end
            end
        end
    end
end
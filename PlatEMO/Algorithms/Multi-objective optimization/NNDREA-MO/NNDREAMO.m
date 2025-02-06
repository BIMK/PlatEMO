classdef NNDREAMO < ALGORITHM
% <2024> <multi> <binary> <large/none> <constrained/none> <sparse/none>
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
                case {'MOKP','Sparse_KP'}
                    Type     = 1;
                    Instance = [Problem.P; Problem.W].';
                case 'Sparse_IS'
                    Type     = 1;
                    Instance = sparse(Problem.Data);
                case 'Sparse_CD'
                    Type    = 2;
                    AdjMat  = Problem.Adj;
                case 'Sparse_CN'
                    Type    = 2;
                    AdjMat  = Problem.A;
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
            % Set the neural network structure
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
            % Get lower and upper bounds of the search space
            lower = repmat(lower,1,Dim);
            upper = repmat(upper,1,Dim);

            %% Generate random population
            % Initialize population weights
            PopWeight = unifrnd(repmat(lower,Problem.N,1),repmat(upper,Problem.N,1));
            % Fully connected neural network for forward propagation
            Output = FCNForward(PopWeight, Instance, s_list);
            % Evaluate the population
            Population = Problem.Evaluation(Output);
            [~,~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,PopWeight,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Problem.FE <= Problem.maxFE * delta
                    % First stage
                    MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                    OffWeight  = OperatorReal(PopWeight(MatingPool, :),lower, upper);
                    Output     = FCNForward(OffWeight, Instance, s_list);
                    Offspring  = Problem.Evaluation(Output);
                    [Population,PopWeight,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[PopWeight;OffWeight], Problem.N);
                else
                    % Second stage
                    for i = 1 : Problem.N
                        drawnow('limitrate');
                        Offspring = OperatorGAhalf(Problem,Population(randperm(end,2)));
                        if all(Offspring.con<=0)
                            [Population,FrontNo] = Reduce([Population,Offspring],FrontNo);
                        end
                    end
                end
            end
        end
    end
end
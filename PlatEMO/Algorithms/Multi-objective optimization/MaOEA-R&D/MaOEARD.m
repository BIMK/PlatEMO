classdef MaOEARD < ALGORITHM
% <many> <real/integer/label/binary/permutation>
% Many-objective evolutionary algorithm based on objective space reduction
% and diversity improvement

%------------------------------- Reference --------------------------------
% Z. He and G. G. Yen, Many-objective evolutionary algorithm: Objective
% space reduction and diversity improvement, IEEE Transactions on
% Evolutionary Computation, 2016, 20(1): 145-160.
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
            W = zeros(Problem.M) + 1e-6;
            W(logical(eye(Problem.M))) = 1;

            %% Generate random population
            Problem.N  = max(ceil(Problem.N/Problem.M)*Problem.M,2*Problem.M);
            Population = Problem.Initialization();

            %% Objective space reduction
            % Reducing the objective space
            while Algorithm.NotTerminated(Population) && Problem.FE < Problem.maxFE/2
                % Classification
                [Subpopulation,Z] = Classification(Population,W);
                % Evolve each subpopulation
                for i = 1 : Problem.M
                    MatingPool = randi(Problem.N/Problem.M,1,Problem.N/Problem.M);
                    Offspring  = OperatorGA(Problem,Subpopulation{i}(MatingPool));
                    Subpopulation{i} = [Subpopulation{i},Offspring];
                    ASF = max((Subpopulation{i}.objs-repmat(Z,length(Subpopulation{i}),1))./repmat(W(i,:),length(Subpopulation{i}),1),[],2);
                    [~,rank] = sort(ASF);
                    Subpopulation{i} = Subpopulation{i}(rank(1:Problem.N/Problem.M));
                end
                Population = [Subpopulation{:}];
            end
            % Identify the target points
            TP = UpdateTP(Population,W);

            %% Diversity improvement
            % Generate random population around TPs
            PopDec     = repmat(TP.decs,Problem.N/Problem.M-1,1);
            PopDec     = PopDec.*(1+randn(size(PopDec))/5);
            Population = [TP,Problem.Evaluation(PopDec)];
            % Evolve
            while Algorithm.NotTerminated(Population)
                MatingPool = randi(Problem.N,1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,TP.objs);
                TP         = UpdateTP([Population,TP],W);
            end
        end
    end
end
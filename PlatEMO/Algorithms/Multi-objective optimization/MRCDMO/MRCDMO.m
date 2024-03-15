classdef MRCDMO < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% Multiregional Co-evolutionary Algorithm


% This function is written by Guopeng Chen

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Changed(Problem,Population)
                    % React to the change
                    [Population] = Reinitialization(Problem,Population,LastPopulation);
                else
                    LastPopulation=Population; %保存上一代的种群
                end
                % Generate offspring randomly
                MatingPool = randperm(Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));

                % Elitism strategy
                UniPop = [Population,Offspring];
                PopObj = UniPop.objs;
                [FrontNo,MaxFNo] = NDSort(PopObj,Problem.N);

                % The number of individuals to be selected in the last
                % non-dominated front
                K = Problem.N - sum(FrontNo<MaxFNo);

                if K ~= 0
                    % Normalization
                    pareto_population = find(FrontNo<MaxFNo);
                    last_population   = find(FrontNo == MaxFNo);
                    Zmin = min(PopObj(FrontNo == MaxFNo,:));
                    Zmax = max(PopObj(FrontNo == MaxFNo,:));
                    S = sum(FrontNo == MaxFNo);
                    MaxFnorm = (PopObj(last_population,:)-repmat(Zmin,S,1))./repmat(Zmax-Zmin,S,1);

                    % Clustering-based reference points generation
                    [Ref] = Reference_Generation( MaxFnorm,Problem.M,K);

                    % Clustering-based environmental selection
                    [reference_population] = Reference_Point_Selection(MaxFnorm,last_population,Ref,K,Problem.M);
                else
                    pareto_population    = find(FrontNo<=MaxFNo);
                    reference_population = [];
                end
                Population = UniPop([pareto_population,reference_population]);
            end
        end
    end
end
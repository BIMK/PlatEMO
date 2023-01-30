classdef CAMOEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% Clustering based adaptive multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% Y. Hua, Y. Jin, K. Hao, A clustering-based adaptive evolutionary
% algorithm for multiobjective optimization with irregular Pareto fronts,
% IEEE Transactions on Cybernetics, 2019, 49(7): 2758-2770.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yicun Hua

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
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
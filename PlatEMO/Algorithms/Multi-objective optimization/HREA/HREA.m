classdef HREA < ALGORITHM
% <multi> <real/integer> <multimodal>
% Hierarchy ranking based evolutionary algorithm
% eps --- 0.3 --- Parameter for quality of the local Pareto front

%------------------------------- Reference --------------------------------
% W. Li, X. Yao, T. Zhang, R. Wang, and L. Wang, Hierarchy ranking method
% for multimodal multi-objective optimization with local Pareto fronts,
% IEEE Transactions on Evolutionary Computation, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenhua Li

    methods
        function main(Algorithm, Problem)
            eps = Algorithm.ParameterSet(0.3);
            p   = 0.5;
            %% Generate random population
            Population          = Problem.Initialization();
            [~,CrowdDis1]       = EnvironmentalSelection(Population,Problem.N);
            [Archive,CrowdDis2] = ArchiveUpdate(Population,Problem.N,eps,0);

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                if Problem.FE >= Problem.maxFE * 0.5 && rand < p
                    MatingPool2 = TournamentSelection(2,round(Problem.N),-CrowdDis2);
                    Offspring   = OperatorGA(Problem,Archive(MatingPool2));
                else
                    MatingPool1 = TournamentSelection(2,round(Problem.N),-CrowdDis1);
                    Offspring   = OperatorGA(Problem,Population(MatingPool1));
                end
                [Population,CrowdDis1] = EnvironmentalSelection([Population,Offspring],Problem.N);
                [Archive,CrowdDis2]    = ArchiveUpdate([Archive,Offspring],Problem.N,eps,Problem.FE/Problem.maxFE);
            end
        end
    end
end
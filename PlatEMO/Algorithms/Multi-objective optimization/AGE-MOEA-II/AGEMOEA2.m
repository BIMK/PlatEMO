classdef AGEMOEA2 < ALGORITHM
% <multi/many> <real/binary/permutation> <constrained/none>

% AGE-MOEA-II

%------------------------------- Reference --------------------------------
% A. Panichella, An Improved Pareto Front Modeling Algorithm for Large-
% scale Many-Objective Optimization, Genetic and Evolutionary Computation 
% Conference (GECCO), 2022, July 9-13, Boston, MA, USA
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
              MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
              Offspring  = OperatorGA(Population(MatingPool));
              [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
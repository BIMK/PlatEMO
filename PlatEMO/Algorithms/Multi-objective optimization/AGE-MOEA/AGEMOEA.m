classdef AGEMOEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Adaptive geometry estimation-based many-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% A. Panichella, An adaptive evolutionary algorithm based on non-euclidean
% geometry for many-objective optimization, Proceedings of the Genetic and
% Evolutionary Computation Conference, 2019, 595-603.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Annibale Panichella

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
              MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
              Offspring  = OperatorGA(Problem,Population(MatingPool));
              [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
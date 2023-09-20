classdef CMME < ALGORITHM
% <many> <real/integer/label/binary/permutation> <constrained>
% Constrained many-objective evolutionary algorithm with enhanced mating and environmental selections

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, L. Wang, and L. Gao, A constrained many-objective
% optimization evolutionary algorithm with enhanced mating and
% environmental selections, IEEE Transactions on Cybernetics, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    methods
        function main(Algorithm,Problem)
            %% Generate Uniform Reference Points
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);

            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                density    = DensityEstimate(Population.objs,W);
                Fitness    = CalFitness2(Population.objs,Population.cons);
                MatingPool = TournamentSelection(2,Problem.N,Fitness,density);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = EnvironmentalSelection([Population,Offspring],W,Problem.N);
            end
        end
    end
end
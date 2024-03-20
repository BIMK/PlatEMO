classdef SPEA2SDE < ALGORITHM
% <many> <real/integer/label/binary/permutation>
% SPEA2 with shift-based density estimation

%------------------------------- Reference --------------------------------
% M. Li, S. Yang, and X. Liu, Shift-based density estimation for
% Pareto-based algorithms in many-objective optimization, IEEE Transactions
% on Evolutionary Computation, 2014, 18(3): 348-365.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
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
            Fitness    = CalFitness(Population.objs);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,Fitness);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
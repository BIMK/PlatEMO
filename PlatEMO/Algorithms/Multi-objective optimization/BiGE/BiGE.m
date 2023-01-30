classdef BiGE < ALGORITHM
% <many> <real/integer/label/binary/permutation>
% Bi-goal evolution

%------------------------------- Reference --------------------------------
% M. Li, S. Yang, and X. Liu, Bi-goal evolution for many-objective
% optimization problems, Artificial Intelligence, 2015, 228: 45-65.
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
            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Estimation(Population.objs,1/Problem.N^(1/Problem.M)));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
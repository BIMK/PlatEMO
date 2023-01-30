classdef MOEAIGDNS < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% Multi-objective evolutionary algorithm based on an enhanced IGD

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zhang, R. Cheng, and Y. Jin, A multi-objective evolutionary
% algorithm based on an enhanced inverted generational distance metric,
% Proceedings of the IEEE Congress on Evolutionary Computation, 2016,
% 5222-5229.
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
            %% Generate the sampling points and random population
            Population = Problem.Initialization();
            Archive    = UpdateArchive(Population,5*Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = randi(Problem.N,1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Archive    = UpdateArchive([Archive,Offspring],5*Problem.N);
                Population = EnvironmentalSelection([Population,Offspring],Archive.objs,Problem.N);
            end
        end
    end
end
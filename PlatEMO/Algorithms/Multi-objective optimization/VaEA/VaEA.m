classdef VaEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Vector angle based evolutionary algorithm

%------------------------------- Reference --------------------------------
% Y. Xiang, Y. Zhou, M. Li, and Z. Chen, A vector angle-based evolutionary
% algorithm for unconstrained many-objective optimization, IEEE
% Transactions on Evolutionary Computation, 2017, 21(1): 131-152.
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
                MatingPool = randi(Problem.N,1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));    
                Population = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
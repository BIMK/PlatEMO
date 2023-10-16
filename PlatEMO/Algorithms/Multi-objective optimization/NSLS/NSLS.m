classdef NSLS < ALGORITHM
% <multi> <real/integer>
% Multiobjective optimization framework based on nondominated sorting and
% local search

%------------------------------- Reference --------------------------------
% B. Chen, W. Zeng, Y. Lin, and D. Zhang, A new local search-based
% multiobjective optimization algorithm, IEEE Transactions on Evolutionary
% Computation, 2015, 19(1): 50-73.
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
                Offspring  = Operator(Problem,Population);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
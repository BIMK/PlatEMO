classdef PESAII < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% Pareto envelope-based selection algorithm II
% div --- 10 --- The number of divisions in each objective

%------------------------------- Reference --------------------------------
% D. W. Corne, N. R. Jerram, J. D. Knowles, and M. J. Oates, PESA-II:
% Region-based selection in evolutionary multiobjective optimization,
% Proceedings of the Annual Conference on Genetic and Evolutionary
% Computation, 2001, 283-290.
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
            %% Parameter setting
            div = Algorithm.ParameterSet(10);

            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs,Problem.N,div);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,div);
            end
        end
    end
end
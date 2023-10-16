classdef GrEA < ALGORITHM
% <many> <real/integer/label/binary/permutation>
% Grid-based evolutionary algorithm
% div --- --- The number of divisions in each objective

%------------------------------- Reference --------------------------------
% S. Yang, M. Li, X. Liu, and J. Zheng, A grid-based evolutionary algorithm
% for many-objective optimization, IEEE Transactions on Evolutionary
% Computation, 2013, 17(5): 721-736.
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
            Div = [0 45 15 10 9 9 8 8 10 12];
            div = Algorithm.ParameterSet(Div(min(Problem.M,10)));

            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs,div);
                Offspring  = OperatorGA(Problem,Population(MatingPool));    
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,div);
            end
        end
    end
end
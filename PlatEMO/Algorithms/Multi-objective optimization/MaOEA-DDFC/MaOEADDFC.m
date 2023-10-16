classdef MaOEADDFC < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Many-objective evolutionary algorithm based on directional diversity and
% favorable convergence
% K --- 5 --- The number of neighbors for estimating density
% L --- 3 --- The number of candidates for convergence-based selection

%------------------------------- Reference --------------------------------
% J. Cheng, G. G. Yen, and G. Zhang, A many-objective evolutionary
% algorithm with enhanced mating and environmental selections, IEEE
% Transactions on Evolutionary Computation, 2015, 19(4): 592-605.
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
            [K,L] = Algorithm.ParameterSet(5,3);

            %% Generate random population
            Population = Problem.Initialization();
            Zmin       = min(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs,Zmin);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Zmin       = min([Zmin;Offspring.objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Zmin,Problem.N,K,L);
            end
        end
    end
end
classdef MOPSO < ALGORITHM
% <multi> <real/integer>
% Multi-objective particle swarm optimization
% div --- 10 --- The number of divisions in each objective

%------------------------------- Reference --------------------------------
% C. A. Coello Coello and M. S. Lechuga, MOPSO: A proposal for multiple
% objective particle swarm optimization, Proceedings of the IEEE Congress
% on Evolutionary Computation, 2002, 1051-1056.
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
            Archive    = UpdateArchive(Population,Problem.N,div);
            Pbest      = Population;

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                REP        = REPSelection(Archive.objs,Problem.N,div);
                Population = OperatorPSO(Problem,Population,Pbest,Archive(REP));
                Archive    = UpdateArchive([Archive,Population],Problem.N,div);
                Pbest      = UpdatePbest(Pbest,Population);
            end
        end
    end
end
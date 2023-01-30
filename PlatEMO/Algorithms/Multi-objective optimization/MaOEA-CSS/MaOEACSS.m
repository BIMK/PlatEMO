classdef MaOEACSS < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Many-objective evolutionary algorithms based on coordinated selection
% strategy
% t --- 0 --- Threshold value in environmental selection

%------------------------------- Reference --------------------------------
% Z. He and G. G. Yen, Many-objective evolutionary algorithms based on
% coordinated selection strategy, IEEE Transactions on Evolutionary
% Computation, 2017, 21(2): 220-233.
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
            t = Algorithm.ParameterSet(0);

            %% Generate random population
            Population = Problem.Initialization();
            Zmin       = min(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs,Zmin);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Zmin       = min([Zmin;Offspring.objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Zmin,t,Problem.N);
            end
        end
    end
end
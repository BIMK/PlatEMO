classdef PICEAg < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Preference-inspired coevolutionary algorithm with goals
% NGoal --- --- Number of goals

%------------------------------- Reference --------------------------------
% R. Wang, R. C. Purshouse, and P. J. Fleming, Preference-inspired
% coevolutionary algorithms for many-objective optimization, IEEE
% Transactions on Evolutionary Computation, 2013, 17(4): 474-494.
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
            NGoal = Algorithm.ParameterSet(100*Problem.M);

            %% Generate random population and goals
            Population = Problem.Initialization();
            Goal       = GeneGoal(Population.objs,NGoal);
            Archive    = UpdateArchive(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                MatingPool = randi(Problem.N,1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Archive    = UpdateArchive([Archive,Offspring],Problem.N);
                newGoal    = GeneGoal([Population.objs;Offspring.objs],NGoal);
                [Population,Goal] = EnvironmentSelection([Population,Offspring],[Goal;newGoal],Problem.N);
            end
        end
    end
end
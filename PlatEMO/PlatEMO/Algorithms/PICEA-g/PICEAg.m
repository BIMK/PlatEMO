function PICEAg(Global)
% <algorithm> <P>
% Preference-inspired coevolutionary algorithm with goals
% NGoal --- --- Number of goals

%------------------------------- Reference --------------------------------
% R. Wang, R. C. Purshouse, and P. J. Fleming, Preference-inspired
% coevolutionary algorithms for many-objective optimization, IEEE
% Transactions on Evolutionary Computation, 2013, 17(4): 474-494.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    NGoal = Global.ParameterSet(100*Global.M);

    %% Generate random population and goals
    Population = Global.Initialization();
    Goal       = GeneGoal(Population.objs,NGoal);
    Archive    = UpdateArchive(Population,Global.N);
    
    %% Optimization
    while Global.NotTermination(Archive)
        MatingPool = randi(Global.N,1,Global.N);
        Offspring  = GA(Population(MatingPool));
        Archive    = UpdateArchive([Archive,Offspring],Global.N);
        newGoal    = GeneGoal([Population.objs;Offspring.objs],NGoal);
        [Population,Goal] = EnvironmentSelection([Population,Offspring],[Goal;newGoal],Global.N);
    end
end
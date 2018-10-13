function PICEAg(Global)
% <algorithm> <O-Z>
% Preference-Inspired Coevolutionary Algorithms for Many-Objective
% Optimization
% NGoal --- --- Number of goals

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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
        Offspring  = Global.Variation(Population(MatingPool));
        Archive    = UpdateArchive([Archive,Offspring],Global.N);
        newGoal    = GeneGoal([Population.objs;Offspring.objs],NGoal);
        [Population,Goal] = EnvironmentSelection([Population,Offspring],[Goal;newGoal],Global.N);
    end
end
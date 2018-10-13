function AGEII(Global)
% <algorithm> <A-G>
% A Fast Approximation-Guided Evolutionary Multi-Objective Algorithm
% epsilon --- 0.1 --- The parameter in grid location calculation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    epsilon = Global.ParameterSet(0.1);

    %% Generate the sampling points and random population
    Population = Global.Initialization();
    Archive    = UpdateArchive(Population,epsilon);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs);
        Offspring  = Global.Variation(Population(MatingPool));
        Archive    = UpdateArchive([Archive,Offspring],epsilon);
        Population = EnvironmentalSelection([Population,Offspring],Archive.objs,Global.N);
    end
end
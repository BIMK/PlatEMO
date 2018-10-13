function MOEAIGDNS(Global)
% <algorithm> <H-N>
% A Multi-objective Evolutionary Algorithm Based on an Enhanced Inverted
% Generational Distance Metric

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the sampling points and random population
    Population = Global.Initialization();
    Archive    = UpdateArchive(Population,5*Global.N);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = randi(Global.N,1,Global.N);
        Offspring  = Global.Variation(Population(MatingPool));
        Archive    = UpdateArchive([Archive,Offspring],5*Global.N);
        Population = EnvironmentalSelection([Population,Offspring],Archive.objs,Global.N);
    end
end
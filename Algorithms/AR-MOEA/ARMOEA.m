function ARMOEA(Global)
% <algorithm> <A-G>
% An Indicator Based Multi-Objective Evolutionary Algorithm with Reference
% Point Adaptation for Better Versatility

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
    W          = UniformPoint(Global.N,Global.M);
    [Archive,RefPoint,Range] = UpdateRefPoint(Population(all(Population.cons<=0,2)).objs,W,[]);
    
    %% Start the interations
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population,RefPoint,Range);
        Offspring  = Global.Variation(Population(MatingPool));
        [Archive,RefPoint,Range] = UpdateRefPoint([Archive;Offspring(all(Offspring.cons<=0,2)).objs],W,Range);
        [Population,Range]       = EnvironmentalSelection([Population,Offspring],RefPoint,Range,Global.N);
    end
end
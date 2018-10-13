function tDEA(Global)
% <algorithm> <O-Z>
% A New Dominance Relation Based Evolutionary Algorithm for Many-Objective
% Optimization

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference points and random population
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    Population   = Global.Initialization();
    [z,znad]     = deal(min(Population.objs),max(Population.objs));

    %% Optimization
    while Global.NotTermination(Population) 
        MatingPool = randi(Global.N,1,Global.N);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,z,znad] = EnvironmentalSelection([Population,Offspring],W,Global.N,z,znad);
    end
end
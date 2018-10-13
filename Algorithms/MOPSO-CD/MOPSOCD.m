function MOPSOCD(Global)
% <algorithm> <H-N>
% An Effective Use of Crowding Distance in Multiobjective Particle Swarm
% Optimization
% operator --- MOPSOCD_operator

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Generate random population
    Population = Global.Initialization();
    Archive    = UpdateArchive(Population,Global.N);
    Pbest      = Population;
    
    %% Optimization
    while Global.NotTermination(Archive)
        Gbest      = Archive(randi(ceil(length(Archive)*0.1),1,Global.N));
        Population = Global.Variation([Population,Pbest,Gbest],Global.N,@MOPSOCD_operator);
        Archive    = UpdateArchive([Archive,Population],Global.N);
        Pbest      = UpdatePbest(Pbest,Population);
    end
end
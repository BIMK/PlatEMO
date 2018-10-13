function MOPSO(Global)
% <algorithm> <H-N>
% Multi-objective particle swarm optimization
% div --- 10 --- The number of divisions in each objective
% operator   --- PSO

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    div = Global.ParameterSet(10);
    
    %% Generate random population
    Population = Global.Initialization();
    Archive    = UpdateArchive(Population,Global.N,div);
    Pbest      = Population;
    
    %% Optimization
    while Global.NotTermination(Archive)
        REP        = REPSelection(Archive.objs,Global.N,div);
        Population = Global.Variation([Population,Pbest,Archive(REP)],Global.N,@PSO);
        Archive    = UpdateArchive([Archive,Population],Global.N,div);
        Pbest      = UpdatePbest(Pbest,Population);
    end
end
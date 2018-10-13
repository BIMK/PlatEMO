function NMPSO(Global)
% <algorithm> <H-N>
% Particle Swarm Optimization with A Balanceable Fitness Estimation for
% Many-objective Optimization Problems
% operator --- NMPSO_operator

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
    Pbest      = Population;
    Archive    = UpdateArchive(Population(NDSort(Population.objs,1)==1),[],Global.N);

    %% Optimization
    while Global.NotTermination(Archive)
        Population = Global.Variation([Population,Pbest,Archive(randi(ceil(length(Archive)/10),1,Global.N))],Global.N,@NMPSO_operator);
        Pbest      = UpdatePbest(Pbest,Population);
        Archive    = UpdateArchive(Archive,Population,Global.N);
        S          = Global.Variation(Archive([1:length(Archive),randi(ceil(length(Archive)/2),1,length(Archive))]),length(Archive),@EAreal);
        Archive    = UpdateArchive(Archive,S,Global.N);
    end
end
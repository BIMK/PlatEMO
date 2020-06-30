function MOPSOCD(Global)
% <algorithm> <M>
% MOPSO with crowding distance

%------------------------------- Reference --------------------------------
% C. R. Raquel and P. C. Naval Jr, An effective use of crowding distance in
% multiobjective particle swarm optimization, Proceedings of the 7th Annual
% Conference on Genetic and Evolutionary Computation, 2005, 257-264.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Generate random population
    Population = Global.Initialization();
    Archive    = UpdateArchive(Population,Global.N);
    Pbest      = Population;
    
    %% Optimization
    while Global.NotTermination(Archive)
        Gbest      = Archive(randi(ceil(length(Archive)*0.1),1,Global.N));
        Population = Operator(Population,Pbest,Gbest);
        Archive    = UpdateArchive([Archive,Population],Global.N);
        Pbest      = UpdatePbest(Pbest,Population);
    end
end
function NMPSO(Global)
% <algorithm> <N>
% Novel multi-objective particle swarm optimization

%------------------------------- Reference --------------------------------
% Q. Lin, S. Liu, Q. Zhu, C. Tang, R. Song, J. Chen, C. A. Coello Coello,
% K. Wong, and J. Zhang, Particle swarm optimization with a balanceable
% fitness estimation for many-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2018, 22(1): 32-46.
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
    Pbest      = Population;
    Archive    = UpdateArchive(Population(NDSort(Population.objs,1)==1),[],Global.N);

    %% Optimization
    while Global.NotTermination(Archive)
        Population = Operator(Population,Pbest,Archive(randi(ceil(length(Archive)/10),1,Global.N)));
        Pbest      = UpdatePbest(Pbest,Population);
        Archive    = UpdateArchive(Archive,Population,Global.N);
        S          = GAhalf(Archive([1:length(Archive),randi(ceil(length(Archive)/2),1,length(Archive))]));
        Archive    = UpdateArchive(Archive,S,Global.N);
    end
end
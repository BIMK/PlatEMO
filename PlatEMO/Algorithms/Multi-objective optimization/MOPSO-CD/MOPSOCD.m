classdef MOPSOCD < ALGORITHM
% <multi> <real/integer>
% MOPSO with crowding distance

%------------------------------- Reference --------------------------------
% C. R. Raquel and P. C. Naval Jr, An effective use of crowding distance in
% multiobjective particle swarm optimization, Proceedings of the Annual
% Conference on Genetic and Evolutionary Computation, 2005, 257-264.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            Archive    = UpdateArchive(Population,Problem.N);
            Pbest      = Population;

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                Gbest      = Archive(randi(ceil(length(Archive)*0.1),1,Problem.N));
                Population = Operator(Problem,Population,Pbest,Gbest);
                Archive    = UpdateArchive([Archive,Population],Problem.N);
                Pbest      = UpdatePbest(Pbest,Population);
            end
        end
    end
end
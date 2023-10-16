classdef NMPSO < ALGORITHM
% <multi/many> <real/integer>
% Novel multi-objective particle swarm optimization

%------------------------------- Reference --------------------------------
% Q. Lin, S. Liu, Q. Zhu, C. Tang, R. Song, J. Chen, C. A. Coello Coello,
% K. Wong, and J. Zhang, Particle swarm optimization with a balanceable
% fitness estimation for many-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2018, 22(1): 32-46.
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
            Pbest      = Population;
            Archive    = UpdateArchive(Population(NDSort(Population.objs,1)==1),[],Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                Population = Operator(Problem,Population,Pbest,Archive(randi(ceil(length(Archive)/10),1,Problem.N)));
                Pbest      = UpdatePbest(Pbest,Population);
                Archive    = UpdateArchive(Archive,Population,Problem.N);
                S          = OperatorGAhalf(Problem,Archive([1:length(Archive),randi(ceil(length(Archive)/2),1,length(Archive))]));
                Archive    = UpdateArchive(Archive,S,Problem.N);
            end
        end
    end
end
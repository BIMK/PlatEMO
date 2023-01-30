classdef MMOPSO < ALGORITHM
% <multi> <real/integer>
% MOPSO with multiple search strategies

%------------------------------- Reference --------------------------------
% Q. Lin, J. Li, Z. Du, J. Chen, and Z. Ming, A novel multi-objective
% particle swarm optimization with multiple search strategies, European
% Journal of Operational Research, 2015, 247(3): 732-744.
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
            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            W = W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));

            %% Generate random population
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);
            Population = Classification(Problem,Population,W,Z);
            Archive    = UpdateArchive(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                [Pbest,Gbest] = GetBest(Archive,W,Z);
                Population    = Operator(Problem,Population,Pbest,Gbest);
                Z             = min([Z;Population.objs],[],1);
                Archive       = UpdateArchive([Archive,Population],Problem.N); 
                S             = OperatorGAhalf(Problem,Archive([1:length(Archive),randi(ceil(length(Archive)/2),1,length(Archive))]));
                Z             = min([Z;S.objs],[],1);
                Archive       = UpdateArchive([Archive,S],Problem.N);
            end
        end
    end
end
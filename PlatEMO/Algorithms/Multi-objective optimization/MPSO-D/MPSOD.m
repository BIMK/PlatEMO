classdef MPSOD < ALGORITHM
% <multi/many> <real/integer>
% Multi-objective particle swarm optimization algorithm based on
% decomposition

%------------------------------- Reference --------------------------------
% C. Dai, Y. Wang, and M. Ye, A new multi-objective particle swarm
% optimization algorithm based on decomposition, Information Sciences,
% 2015, 325: 541-557.
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
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization(2*Problem.N);
            Z          = min(Population.objs,[],1);
            Population = Classification(Problem,Population,W,Z);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                [Parent,Pbest,Gbest] = MatingSelection(Population.objs,B,W,Z);
                Offspring  = Operator(Problem,Population(Parent),Population(Pbest),Population(Gbest));
                Z          = min([Z;Offspring.objs],[],1);
                Population = Classification(Problem,[Population,Offspring],W,Z);
            end
        end
    end
end
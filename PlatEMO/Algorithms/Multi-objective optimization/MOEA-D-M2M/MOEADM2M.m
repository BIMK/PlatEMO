classdef MOEADM2M < ALGORITHM
% <multi> <real/integer>
% MOEA/D based on MOP to MOP
% K --- 10 --- Number of reference vectors

%------------------------------- Reference --------------------------------
% H. Liu, F. Gu, and Q. Zhang, Decomposition of a multiobjective
% optimization problem into a number of simple multiobjective subproblems,
% IEEE Transactions on Evolutionary Computation, 2014, 18(3): 450-455.
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
            %% Parameter setting
            K = Algorithm.ParameterSet(10);

            %% Generate random population
            [W,K]      = UniformPoint(K,Problem.M);
            Problem.N  = ceil(Problem.N/K)*K;
            S          = Problem.N/K;
            Population = Problem.Initialization();
            Population = Associate(Population,W,S);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPoolLocal      = randi(S,S,K) + repmat(0:S:S*(K-1),S,1);
                MatingPoolGlobal     = randi(Problem.N,1,Problem.N);
                rnd                  = rand(S,K) < 0.7;
                MatingPoolLocal(rnd) = MatingPoolGlobal(rnd);
                Offspring  = Operator(Problem,Population,Population(MatingPoolLocal(:)));
                Population = Associate([Population,Offspring],W,S);
            end
        end
    end
end
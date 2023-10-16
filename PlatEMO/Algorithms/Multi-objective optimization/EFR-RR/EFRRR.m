classdef EFRRR < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Ensemble fitness ranking with a ranking restriction scheme
% K --- 2 --- Number of nearest weight vectors

%------------------------------- Reference --------------------------------
% Y. Yuan, H. Xu, B. Wang, B. Zhang, and X. Yao, Balancing convergence and
% diversity in decomposition-based many-objective optimizers, IEEE
% Transactions on Evolutionary Computation, 2016, 20(2): 180-198.
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
            K = Algorithm.ParameterSet(2);

            %% Generate the reference points and random population
            [W,Problem.N]   = UniformPoint(Problem.N,Problem.M);
            Population      = Problem.Initialization();
            [PopObj,z,znad] = Normalization(Population.objs,min(Population.objs),max(Population.objs));
            RgFrontNO       = MaximumRanking(PopObj,W,K);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,RgFrontNO);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,z,znad] = EnvironmentalSelection([Population,Offspring],W,Problem.N,K,z,znad);
            end
        end
    end
end
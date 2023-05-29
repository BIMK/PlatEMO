classdef TSNSGAII < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Two stage NSGA-II

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, and L. Wang, A two-stage evolutionary algorithm with
% balanced convergence and diversity for many-objective optimization, IEEE
% Transactions on Systems, Man, and Cybernetics: Systems, 2022, 52(10):
% 6222-6234.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    methods
        function main(Algorithm,Problem)
            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);

            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,d2] = EnvironmentalSelection(Population,W,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,d2);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                if Problem.FE<(0.8+0.2*(1-3/Problem.M))*Problem.maxFE
                    [Population,FrontNo,d2] = EnvironmentalSelection([Population,Offspring],W,Problem.N);
                else
                    Population = EnvironmentalSelection1([Population,Offspring],W,Problem.N,Problem.M);
                    FrontNo    = NDSort(Population.objs,inf);
                    d2 = DensityEstimate(Population,W);
                end
            end
        end
    end
end
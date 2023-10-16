classdef RPDNSGAII < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Reference point dominance-based NSGA-II

%------------------------------- Reference --------------------------------
% M. Elarbi, S. Bechikh, A. Gupta, L. B. Said, and Y. S. Ong, A new
% decomposition-based NSGA-II for many-objective optimization, IEEE
% Transactions on Systems, Man, and Cybernetics: Systems, 2018, 48(7):
% 1191-1210.
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
            %% Generate the reference points and random population
            [RPSet,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population        = Problem.Initialization();
            [~,FrontNo,d2]    = EnvironmentalSelection(Population,RPSet,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population) 
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,d2);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,d2] = EnvironmentalSelection([Population,Offspring],RPSet,Problem.N);
            end
        end
    end
end
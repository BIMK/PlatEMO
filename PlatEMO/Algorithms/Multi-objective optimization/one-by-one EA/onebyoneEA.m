classdef onebyoneEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Many-objective evolutionary algorithm using a one-by-one selection
% strategy

%------------------------------- Reference --------------------------------
% Y. Liu, D. Gong, J. Sun, and Y. Jin, A many-objective evolutionary
% algorithm using a one-by-one selection strategy, IEEE Transactions on
% Cybernetics, 2017, 47(9): 2689-2702.
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
            % Ideal point
            zmin = min(Population.objs,[],1);
            % Rank of each solution in one-by-one selection
            Rank = ones(1,Problem.N);
            % Distribution threshold
            zeta = 1;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs,Rank);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,Rank,zeta,zmin] = EnvironmentalSelection([Population,Offspring],zeta,zmin,Problem.N);
            end
        end
    end
end
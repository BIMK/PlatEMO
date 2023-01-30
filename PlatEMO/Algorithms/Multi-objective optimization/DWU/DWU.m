classdef DWU < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% Dominance-weighted uniformity multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% G. Moreira and L. Paquete, Guiding under uniformity measure in the
% decision space, Proceedings of the IEEE Latin American Conference on
% Computational Intelligence, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Gladston Moreira

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [Population,FrontNo] = EnvironmentalSelection(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,] = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
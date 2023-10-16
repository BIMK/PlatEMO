classdef GDE3 < ALGORITHM
% <multi> <real/integer> <constrained/none>
% Generalized differential evolution 3

%------------------------------- Reference --------------------------------
% S. Kukkonen and J. Lampinen, GDE3: The third evolution step of
% generalized differential evolution, Proceedings of the IEEE Congress on
% Evolutionary Computation, 2005, 443-450.
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

            %% Optimization
            while Algorithm.NotTerminated(Population)
                Offspring  = OperatorDE(Problem,Population,Population(randi(Problem.N,1,Problem.N)),Population(randi(Problem.N,1,Problem.N)));
                Population = EnvironmentalSelection(Population,Offspring,Problem.N);
            end
        end
    end
end
classdef MyODEMR < ALGORITHM
% <multi/many> <real/integer>
% Many-objective differential evolution with mutation restriction
% nP --- 500 --- Number of reference points for IGD calculation

%------------------------------- Reference --------------------------------
% R. Denysiuk, L. Costa, and I. E. Santo, Many-objective optimization using
% differential evolution with variable-wise mutation restriction,
% Proceedings of the Annual Conference on Genetic and Evolutionary
% Computation, 2013, 591-598.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            nP = Algorithm.ParameterSet(500);

            %% Generate hyperplane 
            P = UniformPoint(nP,Problem.M);

            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                Offspring  = Operator(Problem,Population(1:Problem.N),Population(randi(Problem.N,1,Problem.N)),Population(randi(Problem.N,1,Problem.N)));
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,P);
            end
        end
    end
end
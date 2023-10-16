classdef EMyOC < ALGORITHM
% <multi/many> <real/integer>
% Evolutionary many-objective optimization algorithm with clustering-based
% selection

%------------------------------- Reference --------------------------------
% R. Denysiuk, L. Costa, and I. E. Santo, Clustering-based selection for
% evolutionary many-objective optimization, Proceedings of the
% International Conference on Parallel Problem Solving from Nature, 2014,
% 538-547.
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
            %% Generate random population
            % Initialize population
            Population = Problem.Initialization();
            % Initialize ideal point
            Z = min(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                Offspring  = Operator(Problem,Population,Population(randi(Problem.N,1,Problem.N)),Population(randi(Problem.N,1,Problem.N)));
                Z          = min(Z,min(Offspring.objs,[],1));
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z);
            end
        end
    end
end
classdef PeEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Pareto front shape estimation based evolutionary algorithm

%------------------------------- Reference --------------------------------
% L. Li, G. G. Yen, A. Sahoo, L. Chang, and T. Gu, On the estimation of
% pareto front and dimensional similarity in many-objective evolutionary
% algorithm, Information Sciences, 2021, 563: 375-400.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Li Li

	methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = EnvironmentalSelection([Population,Offspring],Problem.N); 
            end
        end
    end
end
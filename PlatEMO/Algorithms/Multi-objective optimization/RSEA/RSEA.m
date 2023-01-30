classdef RSEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Radial space division based evolutionary algorithm

%------------------------------- Reference --------------------------------
% C. He, Y. Tian, Y. Jin, X. Zhang, and L. Pan, A radial space division
% based evolutionary algorithm for many-objective optimization, Applied
% Soft Computing, 2017, 61: 603-621.
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
            Range      = inf(2,Problem.M);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                Range(1,:) = min([Range(1,:);Population.objs],[],1);
                Range(2,:) = max(Population(NDSort(Population.objs,1)==1).objs,[],1);
                MatingPool = MatingSelection(Population.objs,Range,ceil(Problem.N/2)*2);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = EnvironmentalSelection(Problem,[Population,Offspring],Range,Problem.N);
            end
        end
    end
end
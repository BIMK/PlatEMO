classdef CPSMOEA < ALGORITHM
% <multi> <real/integer> <expensive>
% Classification and Pareto domination based multi-objective evolutionary
% algorithm
% M --- 3 --- Number of generated offsprings for each solution

%------------------------------- Reference --------------------------------
% J. Zhang, A. Zhou, and G. Zhang, A classification and Pareto domination
% based multiobjective evolutionary algorithm, Proceedings of the IEEE
% Congress on Evolutionary Computation, 2015, 2883-2890.
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
            M = Algorithm.ParameterSet(3);

            %% Generate random population
            Population   = Problem.Initialization();
            [Pgood,Pbad] = NDS(Population,floor(Problem.N/2));

            %% Optimization
            while Algorithm.NotTerminated(Population)
                KNN(Pgood.decs,Pbad.decs);
                Offspring  = Operator(Problem,Population,M);
                Population = NDS([Population,Offspring],Problem.N);
                FrontNo    = NDSort(Offspring.objs,1);
                Pgood      = NDS([Pgood,Offspring(FrontNo==1)],floor(Problem.N/2));
                Pbad       = NDS([Pbad,Offspring(FrontNo~=1)],floor(Problem.N/2));
            end
        end
    end
end
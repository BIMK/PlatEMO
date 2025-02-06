classdef WASFGA < ALGORITHM
% <2015> <multi> <real/integer/label/binary/permutation>
% Weighting achievement scalarizing function genetic algorithm
% Point --- --- Preferred point

%------------------------------- Reference --------------------------------
% A. B. Ruiz, R. Saborido, and M. Luque. A preference-based evolutionary
% algorithm for multiobjective optimization: the weighting achievement
% scalarizing function genetic algorithm. Journal of Global Optimization,
% 2015, 62: 101-129.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            Point = Algorithm.ParameterSet(zeros(1,Problem.M)+0.5,0.00001);
            ro    = 0.0001;

            %% Generate random population
            Population = Problem.Initialization();

            %% Generate a sample of weight vectors
            [n,~] = size(Population.objs);
            if Problem.M == 2
                Vectors = generateWeightVectors2(n, 0.001);
            else
                [Vectors,Problem.N] = UniformPoint(Problem.N,Problem.M);
            end
            FrontNo  = WASFGASort(Vectors, Population.objs, inf, Point,ro);
            CrowdDis = CrowdingDistance(Population.objs,FrontNo);
            [v,~] = size(Vectors);
            if v >= n
                nsort = 2;
            else
                nsort = floor(n/v) + 1;
            end

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelectionW(Vectors, [Population,Offspring],nsort,Point,ro);
            end
        end
    end
end
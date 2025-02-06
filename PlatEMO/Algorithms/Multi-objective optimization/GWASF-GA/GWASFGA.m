classdef GWASFGA < ALGORITHM
% <2017> <multi> <real/integer/label/binary/permutation>
% Global weighting achievement scalarizing function genetic algorithm

%------------------------------- Reference --------------------------------
% R. Saborido, A. B. Ruiz, and M. Luque. Global WASF-GA: An evolutionary
% algorithm in multiobjective optimization to approximate the whole Pareto
% optimal front. Evolutionary computation, 2017, 25(2): 309-349.
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
            ro  = 0.0001;
            eps = 0.01;
            
            %% Generate random population
            Population = Problem.Initialization();

            %% Generate a sample of weight vectors
            [n,col] = size(Population.objs);
            if Problem.M == 2
                Vectors = generateWeightVectors2(n, 0.001);
            else
                [Vectors,Problem.N] = UniformPoint(Problem.N,Problem.M );
            end
            [v,~] = size(Vectors);
            if v >= n
                nsort = 2;
            else
                nsort = floor(n/v) + 1;
            end
            nadir = zeros(1, col);
            Utop  = zeros(1, col);
            disp(Population.objs);
            disp(size(Population.objs))
            A = Population.objs;

            %% Initialize Nadir and Utopian points.
            for i = 1:col
                Maxs = max(A(:, i)) + eps;
                Mins = min(A(:, i)) - eps;
                nadir(i) = Maxs;
                Utop(i)  = Mins;
            end
            FrontNo  = GWASFGASort(Vectors, Population.objs,Utop,nadir, nsort,ro, eps);
            CrowdDis = CrowdingDistance(Population.objs,FrontNo);
                                
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelectionGW(Vectors, [Population,Offspring], Utop,nadir, nsort,ro, eps);
                P = Population.objs;
                % Check if nadir and Utopian points have changed after each
                % generation of the algorithm
                for i = 1 : col
                    Maxs = max(P(:,i)) + eps;
                    Mins = min(P(:,i)) - eps;
                    if Utop(i) > Mins
                        Utop(i) = Mins;
                    end
                    if nadir(i) < Maxs        
                        nadir(i) = Maxs;
                    end
                end
            end
        end
    end
end
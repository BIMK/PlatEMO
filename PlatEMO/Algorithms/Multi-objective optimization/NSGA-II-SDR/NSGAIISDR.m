classdef NSGAIISDR < ALGORITHM
% <many> <real/integer/label/binary/permutation>
% NSGA-II with strengthened dominance relation

%------------------------------- Reference --------------------------------
% Y. Tian, R. Cheng, X. Zhang, Y. Su, and Y. Jin, A strengthened dominance
% relation considering convergence and diversity for evolutionary many-
% objective optimization, IEEE Transactions on Evolutionary Computation,
% 2019, 23(2): 331-345.
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
            zmin       = min(Population.objs,[],1);
            zmax       = max(Population.objs,[],1);
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N,zmin,zmax);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                zmin       = min([zmin;Offspring.objs],[],1);
                zmax       = max(Population(FrontNo==1).objs,[],1);
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N,zmin,zmax);
            end
        end
    end
end
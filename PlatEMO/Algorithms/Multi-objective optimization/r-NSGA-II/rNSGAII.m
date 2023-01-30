classdef rNSGAII < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% r-dominance based NSGA-II
% Points ---     --- Set of preferred points
% W      ---     --- Set of weight vectors for preferred points
% delta  --- 0.1 --- Non-r-dominance threshold

%------------------------------- Reference --------------------------------
% L. B. Said, S. Bechikh, and K. Ghedira, The r-dominance: A new dominance
% relation for interactive evolutionary multicriteria decision making, IEEE
% Transactions on Evolutionary Computation, 2010, 14(5): 801-818.
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
            [Points,W,delta] = Algorithm.ParameterSet(zeros(1,Problem.M)+0.5,ones(1,Problem.M),0.1);

            %% Generate random population
            Population = Problem.Initialization();
            FrontNo    = NrDSort(Population.objs,inf,Points,W,1-(1-delta)*Problem.FE/Problem.maxFE);
            CrowdDis   = CrowdingDistance(Population.objs,FrontNo);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N,Points,W,1-(1-delta)*Problem.FE/Problem.maxFE);
            end
        end
    end
end
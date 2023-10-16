classdef DNNSGAII < ALGORITHM
% <multi> <real/integer> <multimodal>
% Decision space based niching NSGA-II

%------------------------------- Reference --------------------------------
% J. Liang, C. Yue, and B. Qu. Multimodal multi-objective optimization: A
% preliminary study, Proceedings of the IEEE Congress on Evolutionary
% Computation, 2016, 2454-2461.
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
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection_Mod(round(Problem.N/2),round(Problem.N/2),Population.decs,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
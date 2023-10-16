classdef NSGAIIARSBX < ALGORITHM
% <multi> <real/integer> <constrained/none>
% NSGA-II with adaptive rotation based simulated binary crossover

%------------------------------- Reference --------------------------------
% L. Pan, W. Xu, L. Li, C. He, and R. Cheng, Adaptive simulated binary
% crossover for rotated multi-objective optimization, Swarm and
% Evolutionary Computation, 2021, 60: 100759.
%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            B  = eye(Problem.D);
            m  = 0.5*(Problem.upper - Problem.lower);
            ps = 0.5;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = ARSBX(Problem,Population(MatingPool),{B,m,ps});
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
                [B,m,ps,Population] = UpdateParameter(Problem,Population);
            end
        end
    end
end
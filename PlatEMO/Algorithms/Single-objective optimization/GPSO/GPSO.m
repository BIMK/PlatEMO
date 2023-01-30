classdef GPSO < ALGORITHM
% <single> <real> <large/none> <constrained/none>
% Gradient based particle swarm optimization algorithm

%------------------------------- Reference --------------------------------
% M. M. Noel, A new gradient based particle swarm optimization algorithm
% for accurate computation of global minimum, Applied Soft Computing, 2012,
% 12: 353-359.
%------------------------------- Copyright --------------------------------
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
            Pbest      = Population;
            [~,best]   = min(FitnessSingle(Pbest));
            Gbest      = LocalSearch(Problem,Pbest(best));
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                Population     = OperatorPSO(Problem,Population,Pbest,Gbest);
                replace        = FitnessSingle(Pbest) > FitnessSingle(Population);
                Pbest(replace) = Population(replace);
                [~,best]       = min(FitnessSingle(Pbest));
                Gbest          = LocalSearch(Problem,Pbest(best));
            end
        end
	end
end
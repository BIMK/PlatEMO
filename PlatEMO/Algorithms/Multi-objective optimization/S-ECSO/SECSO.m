classdef SECSO < ALGORITHM
% <multi> <real> <large/none> <sparse>
% Enhanced competitive swarm optimizer for sparse optimization
    
%------------------------------- Reference --------------------------------
% X. Wang, K. Zhang, J. Wang, and Y. Jin, An enhanced competitive swarm
% optimizer with strongly convex sparse operator for large-scale
% multi-objective optimization, IEEE Transactions on Evolutionary
% Computation, 2022, 26(5): 859-871.
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
            % Used in ECSO
            subswarm       = floor(Problem.N/3);         
            subswarm_index = [1,1+subswarm,Problem.N-subswarm];
            % Used in SCSparse 
            [LMax,LMin] = Algorithm.ParameterSet(0.35,0);
            step        = (LMax-LMin).*(Problem.upper-Problem.lower)/((Problem.maxFE/Problem.N/3)-1);
            lamb        = LMax.*(Problem.upper-Problem.lower);
            %% Population initialization
            x = rand(Problem.N,Problem.D);
            v = rand(Problem.N,Problem.D);
            x = x.*(repmat(Problem.upper,Problem.N,1)-repmat(Problem.lower,Problem.N,1)) + repmat(Problem.lower,Problem.N,1);
            v = v.*(repmat(Problem.upper,Problem.N,1)-repmat(Problem.lower,Problem.N,1)) + repmat(Problem.lower,Problem.N,1);
            Population  = Problem.Evaluation(x);
            Population1 = Population; % 'Population1'is the population(x), 'Population' is A in the paper
            
            %% Optimization
            iter = 0;
            while Algorithm.NotTerminated(Population)
                iter = iter + 1;
                if iter > Problem.maxFE/Problem.N/3-1
                    Problem.FE = Problem.maxFE;
                end
                [Population1]      = SCSparse(Problem,Population1,lamb);
                [Population,gBest] = A_get(Problem,Population1,Population,iter);
                [Population1,v]    = ECSO(Problem,Population1,v,gBest,subswarm_index);
                lamb = lamb - step;     % Update lamb in SCSparse
            end
        end
    end
end
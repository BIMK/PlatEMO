classdef SAMSO < ALGORITHM
% <single> <real/integer> <large/none> <expensive>
% Multiswarm-assisted expensive optimization

%------------------------------- Reference --------------------------------
% F. Li, X. Cai, L. Gao, and W. Shen, A surrogate-assisted multiswarm
% optimization algorithm for high-dimensional computationally expensive
% problems, IEEE Transactions on Cybernetics, 2021, 51(3): 1390-1402.
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
            assert(~isempty(ver('Optim')),'The execution of SAMPSO requires the Optimization Toolbox.');
            
            %% Parameter setting
            [Wnc,Pr] = Algorithm.ParameterSet(1,0.5);
            eta      = min(sqrt(0.001^2*Problem.D),5e-4*min(Problem.upper-Problem.lower));
            
            %% Initialze DB
            if Problem.D > 50
                N = 80;
                K = 2*Problem.D;
            else
                N = 40;
                K = N;
            end
            PopDec = UniformPoint(K,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,K,1).*PopDec+repmat(Problem.lower,K,1));
            
            %% Initialize swarm
            % Determine position
            [~,idx]  = sort(Population.objs,'ascend');
            Select   = idx(1:N);
            Position = [Population(Select).decs,Population(Select).objs];
            % velocity
            Vmax     = 0.5*(Problem.upper-Problem.lower);
            Vmin     = -0.5*Vmax;
            Velocity = rand(N,Problem.D).*(repmat(Vmax-Vmin,N,1)) + repmat(Vmin,N,1);
            % Pbest and Gbest
            Pbest    = Position;
            Gbest    = Position(1,:);
            maxFES   = Problem.maxFE - K;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Build RBF surrogate model
                [model,~] = rbf_build(Population.decs,Population.objs);

                % Find the minimum of the surrogate
                srgtMin = FindOpt(model,Population,Problem.upper,Problem.lower);

                % Calculate distance
                dist  = pdist2(Population.decs,srgtMin);
                dxRBF = min(dist);
                if dxRBF > eta
                    optSrgt    = Problem.Evaluation(srgtMin);
                    Population = [Population,optSrgt];
                    if optSrgt.objs < Gbest(:,end)
                        [model,~] = rbf_build(Population.decs,Population.objs);
                        Gbest     = [optSrgt.decs,optSrgt.objs];
                    end
                end
                currFES = Problem.FE - K;
                [Population,Position,Velocity,Gbest,Pbest] = UpdatePosition(Problem,Population,Position,Velocity,Pbest,Gbest,currFES,maxFES,Wnc,Pr,model,eta);
            end
        end
    end
end
classdef dMOPSO < ALGORITHM
% <multi> <real/integer>
% MOPSO based on decomposition
% Ta --- 2 --- Age threshold

%------------------------------- Reference --------------------------------
% S. Z. Martinez and C. A. Coello Coello, A multi-objective particle swarm
% optimizer based on decomposition, Proceedings of the Annual Conference on
% Genetic and Evolutionary Computation, 2011, 69-76.
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
            Ta = Algorithm.ParameterSet(2);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);

            %% Generate random population
            Population = Problem.Initialization();
            Age        = zeros(1,Problem.N);         	% Age of each particle
            Z          = min(Population.objs,[],1);     % Ideal point
            Pbest      = Population;                    % Personal best
            Gbest      = UpdateGbest(W,Population,Z);   % Global best

            %% Optimization
            while Algorithm.NotTerminated(Gbest)
                % Generate offsprings
                Population(Age<Ta) = Operator(Problem,Population(Age<Ta),Pbest(Age<Ta),Gbest(Age<Ta));
                for i = find(Age>=Ta)
                    temp          = randn(1,Problem.D);
                    Population(i) = Problem.Evaluation((Gbest(i).dec-Pbest(i).dec)/2+temp.*abs(Gbest(i).dec-Pbest(i).dec));
                    Age(i)        = 0;
                end

                % Update Z
                Z = min(Z,min(Population.objs,[],1));

                % Update personal best
                PbeObj    = Pbest.objs - repmat(Z,Problem.N,1);
                PopObj    = Population.objs - repmat(Z,Problem.N,1);
                normW     = sqrt(sum(W.^2,2));
                normPbe   = sqrt(sum(PbeObj.^2,2));
                normPop   = sqrt(sum(PopObj.^2,2));
                CosinePbe = sum(PbeObj.*W,2)./normW./normPbe;
                CosinePop = sum(PopObj.*W,2)./normW./normPop;
                g_old     = normPbe.*CosinePbe + 5*normPbe.*sqrt(1-CosinePbe.^2);
                g_new     = normPop.*CosinePop + 5*normPop.*sqrt(1-CosinePop.^2);
                Pbest(g_new<=g_old) = Population(g_new<=g_old);
                Age(g_new<=g_old)   = 0;
                Age(g_new>g_old)    = Age(g_new>g_old) + 1;

                % Update global best
                Gbest = UpdateGbest(W,[Gbest,Population],Z);
            end
        end
    end
end
classdef MASSRMADE < ALGORITHM
% <2026> <single> <real>
% Multi-agent differential evolution with MASSR states

%------------------------------- Reference --------------------------------
% X. Liu, Y. Tian, S. Yang, Z. Sheng, and X. Zhang. A multi-agent
% self-supervised state representation framework for automated algorithm
% configuration. IEEE Transactions on Evolutionary Computation, 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate the initial population
            Population = Problem.Initialization();
            N          = Problem.N;

            %% Initialize representation learning and DE settings
            warmFE  = ceil(0.2*Problem.maxFE);
            encoder = MASSREncoder(Problem.D+Problem.M,'MADE');
            encoder.collect(MASSRFeatures.observation(Population,Problem), ...
            MASSRFeatures.guidance(Population,[],Problem,0,0,0));
            controller = [];
            stagnation = 0;
            window     = max(1,0.1*Problem.maxFE/N);
            parameters = [0.30,0.95;0.50,0.90;0.80,0.60;1.00,0.30];

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Train MASSR once at the random-collection boundary
                if ~encoder.Ready && Problem.FE == warmFE
                    encoder.fit();
                    controller = MASSRMADEPolicy();
                end

                % Configure mating, the DE operator, and its parameters
                if encoder.Ready
                    state    = encoder.encode(MASSRFeatures.observation(Population,Problem));
                    progress = (Problem.FE-warmFE)/max(Problem.maxFE-warmFE-N,1);
                    action   = controller.act(state,true(4,3),progress);
                    limit    = Problem.maxFE;
                else
                    action = [randi(3);randi(4);randi(4)];
                    limit  = warmFE;
                end

                % Generate offspring using the selected mating strategy
                LastPopulation = Population;
                count          = min(N,limit-Problem.FE);
                order          = randperm(N,count);
                fitness        = double(FitnessSingle(Population));
                Offspring      = Population([]);
                for i = order
                    parents          = MASSRMADEParents(Population,fitness,Problem,i,action(1));
                    Offspring(end+1) = MASSRMADEOperator(Problem,parents,action(2), ...
                        parameters(action(3),1),parameters(action(3),2));
                end

                % Replace improved parents and compute the team reward
                changed             = FitnessSingle(Offspring) < fitness(order);
                indices             = order(changed);
                Population(indices) = Offspring(changed);
                old                 = double(LastPopulation.objs);
                now                 = double(Population.objs);
                replacement         = sum(changed)/N;
                reward              = min(old) - min(now) + mean(old) - mean(now) + replacement;
                improvement         = MASSRFeatures.change(min(old),min(now));
                if improvement > 1e-6
                    stagnation = 0;
                else
                    stagnation = stagnation + 1;
                end

                % Store a policy transition or a representation sample
                X = MASSRFeatures.observation(Population,Problem);
                if encoder.Ready
                    next = encoder.encode(X);
                    controller.remember(state,action,reward,next,Problem.FE>=Problem.maxFE,true(4,3));
                    controller.train();
                else
                    encoder.collect(X,MASSRFeatures.guidance(Population,LastPopulation,Problem, ...
                    improvement,stagnation/window,replacement));
                end
            end
        end
    end
end
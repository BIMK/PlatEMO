classdef ABC < ALGORITHM
% <single> <real/integer> <large/none> <constrained/none>
% Artificial bee colony algorithm
% limit --- 20 --- The number of trials for releasing a food source

%------------------------------- Reference --------------------------------
% D. Karaboga, An idea based on honey bee swarm for numerical optimization,
% Erciyes University, Tech. Rep. tr06, 2005.
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
            limit = Algorithm.ParameterSet(20);
            
            %% Generate random population
            Population = Problem.Initialization();
            Limit      = zeros(1,Problem.N);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Employed bees
                Pdec      = Population.decs;
                Odec      = Pdec + (rand(size(Pdec))*2-1).*(Pdec-Pdec(randi(end,1,end),:));
                Offspring = Problem.Evaluation(Odec);
                replace             = FitnessSingle(Population) > FitnessSingle(Offspring);
                Population(replace) = Offspring(replace);
                Limit(~replace)     = Limit(~replace) + 1;
                
                % Onlooker bees
                Q         = RouletteWheelSelection(Problem.N,exp(Population.objs/mean(abs(Population.objs+1e-6))));
                Pdec      = Population.decs;
                Odec      = Pdec(Q,:) + (rand(size(Pdec))*2-1).*(Pdec(Q,:)-Pdec(randi(end,1,end),:));
                Offspring = Problem.Evaluation(Odec);
                replace             = FitnessSingle(Population) > FitnessSingle(Offspring);
                Population(replace) = Offspring(replace);
                Limit(~replace)     = Limit(~replace) + 1;

                % Scout bees
                Q = Limit > limit;
                if any(Q)
                    Population(Q) = Problem.Initialization(sum(Q));
                    Limit(Q)      = 0;
                end
            end
        end
    end
end
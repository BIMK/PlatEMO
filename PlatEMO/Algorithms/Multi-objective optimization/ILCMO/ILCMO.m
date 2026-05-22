classdef ILCMO < ALGORITHM
% <2026> <multi/many> <real/integer> <large/none> <constrained/none>
% Indicator-based evolutionary algorithm for large-scale constrained multi-objective optimization

%------------------------------- Reference --------------------------------
% X. Zhong, X. Yao, K. Qiao, D. Gong, and Y. Jin. An indicator-based
% evolutionary algorithm for large-scale constrained multiobjective
% optimization. IEEE Transactions on Evolutionary Computation, 2026, 30(1):
% 271-285.
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
            %% Parameter setting
            [typeOfGroups,numberOfGroups] = Algorithm.ParameterSet(2,4);
    
            %% Generate the reference points and random population
            W             = UniformPoint(Problem.N,Problem.M);
            Population{1} = Problem.Initialization();
            Population{2} = Problem.Initialization();
            Fmin          = min(Population{1}.objs,[],1);
            Fitness{1}    = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{2}    = CalFitness(Population{2}.objs);
    
            %% Calculate the initial dynamic constraint boundary
            CPopulation = [Population{1},Population{2}];
            VAR0        = max(sum(max(CPopulation.cons,0),2));
            if VAR0 == 0
                VAR0 = 1;
            end
            X = 0;
    
            %% Optimization
            while Algorithm.NotTerminated([Population{1}])
                % Udate the epsilon value
                cp   = (-log(VAR0)-6)/log(1-0.5);
                VAR1 = VAR0*(1-X)^cp;
                rf   = sum(sum(max(Population{1}.cons,0),2)<1e-6)/length(Population{1});
                VAR2 = VAR1*(1-rf);
                X    = X+1/(Problem.maxFE/(2*Problem.N));
                % Offspring generation
                P = 1-(2/(1+exp(-5*Problem.FE/Problem.maxFE))-1);
                for i = 1 : 2
                    valOffspring{i} = VGDE_main(Problem, Population{i}, Fitness{i},numberOfGroups, typeOfGroups,P);
                end
                % Update ideal points
                Offspring = [valOffspring{1},valOffspring{2}];
                Fmin      = min([Fmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                % Environmental selection
                [Population{1},Fitness{1},~] = Main_task_Update([Offspring,Population{1}],Problem.N,Fmin);
                [Population{2},Fitness{2},~] = DI_Update([Offspring,Population{1},Population{2}],Problem.N,VAR1,VAR2,W);
            end
        end
    end
end
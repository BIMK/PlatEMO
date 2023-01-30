classdef BSPGA < ALGORITHM
% <single> <binary> <large/none> <constrained/none>
% Binary space partition tree based genetic algorithm
% proC   ---  0.5 --- Probability of crossover
% proM   ---    1 --- Expectation of the number of mutated variables
% lambda --- 0.05 --- Tree based learning probability

%------------------------------- Reference --------------------------------
% Y. Su, N. Guo, Y. Tian, and X. Zhang, A non-revisiting genetic algorithm
% based on a novel binary space partition tree, Information Sciences, 2020,
% 512: 661-674.
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
            [proC,proM,lambda] = Algorithm.ParameterSet(0.5,1,0.05);
            
            %% Generate random population
            [Population,T,best] = BSPTreeConstruction(Problem,randi([0,1],Problem.N,Problem.D),[],[]);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                OffDec = OperatorGA(Problem,Population.decs,{proC,0,proM,0});
                OffDec = BSPTreeLearning(OffDec,best,lambda);
                [Offspring,T,best] = BSPTreeConstruction(Problem,OffDec,T,best);
                Population = [Population,Offspring];
                [~,rank]   = sort(FitnessSingle(Population));
                Population = Population(rank(1:Problem.N));
            end
        end
    end
end

classdef DE < ALGORITHM
% <single> <real/integer> <large/none> <constrained/none>
% Differential evolution
% CR --- 0.9 --- Crossover constant
% F  --- 0.5 --- Mutation factor

%------------------------------- Reference --------------------------------
% R. Storn and K. Price, Differential evolution-a simple and efficient
% heuristic for global optimization over continuous spaces, Journal of
% Global Optimization, 1997, 11(4): 341-359.
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
            [CR,F] = Algorithm.ParameterSet(0.9,0.5);
            
            %% Generate random population
            Population = Problem.Initialization();
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,2*Problem.N,FitnessSingle(Population));
                Offspring  = OperatorDE(Problem,Population,Population(MatingPool(1:end/2)),Population(MatingPool(end/2+1:end)),{CR,F,0,0});
                replace             = FitnessSingle(Population) > FitnessSingle(Offspring);
                Population(replace) = Offspring(replace);
            end
        end
    end
end
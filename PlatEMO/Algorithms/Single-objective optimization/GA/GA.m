classdef GA < ALGORITHM
% <single> <real/integer/label/binary/permutation> <large/none> <constrained/none>
% Genetic algorithm
% proC ---  1 --- Probability of crossover
% disC --- 20 --- Distribution index of simulated binary crossover
% proM ---  1 --- Expectation of the number of mutated variables
% disM --- 20 --- Distribution index of polynomial mutation

%------------------------------- Reference --------------------------------
% J. H. Holland, Adaptation in Natural and Artificial Systems, MIT Press,
% 1992.
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
            [proC,disC,proM,disM] = Algorithm.ParameterSet(1,20,1,20);
            
            %% Generate random population
            Population = Problem.Initialization();
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FitnessSingle(Population));
                Offspring  = OperatorGA(Problem,Population(MatingPool),{proC,disC,proM,disM});
                Population = [Population,Offspring];
                [~,rank]   = sort(FitnessSingle(Population));
                Population = Population(rank(1:Problem.N));
            end
        end
    end
end
classdef MONAS < ALGORITHM
% <2025> <multi> <real/integer/label/binary/permutation> <multimodal>
% Multi-objective neural architecture search

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, B. Xue, M. Zhang, and Y. Jin. An evolutionary framework
% for multi-objective neural architecture search. IEEE Transactions on
% Evolutionary Computation, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            
            %% calculate fitness of populations
            [Fitness1,D_Dec] = CalFitness(Population1.objs,Population1.decs);
            [Fitness2,D]     = CalFitnessDecAux(Population2.objs,Population2.decs);

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                if Problem.FE <= Problem.maxFE/2
                    MatingPool1 = TournamentSelection(2,Problem.N,D_Dec,Fitness1);
                    MatingPool2 = TournamentSelection(2,Problem.N,D,Fitness2);
                    Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                    Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                else
                    MatingPool1 = TournamentSelection(2,2*Problem.N,D_Dec,Fitness1);
                    MatingPool2 = TournamentSelection(2,2*Problem.N,D,Fitness2);
                    Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                    Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                end
                Offspring = [Offspring1,Offspring2];
                [Population1,Fitness1,D_Dec] = EnvironmentalSelection([Population1,Offspring],Problem.N);
                [Population2,Fitness2,D]     = EnvironmentalSelectionAux([Population2,Offspring],Problem.N);
            end
        end
    end
end
classdef CMMO < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <multimodal>
% Coevolutionary multi-modal multi-objective optimization framework

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, L. Wang, and L. Gao, Balancing convergence and
% diversity in objective and decision spaces for multimodal multi-objective
% optimization, IEEE Transactions on Emerging Topics in Computational
% Intelligence, 2023, 7(2): 474-486.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
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
            
            %% parameters for epsilon and initialization of epsilon
            G = ceil(Problem.maxFE/(2 * Problem.N));
            gen = 1;
            Tc = 0.8 * G;
            tao = 0.1;
            epsilon_0 = sum(max(Population2.objs,[],1),2);
            threshold = 0.1;
            epsilon_k = epsilon_0;
            
            %% calculate fitness of populations
            [Fitness1,D_Dec,D_Pop] = CalFitness(Population1.objs,Population1.decs);
            [Fitness2,D] = CalFitnessDecEpsilon(Population2.objs,Population2.decs,epsilon_k);

            %% Optimization
            while Algorithm.NotTerminated(Population1)         
                % The value of e(k) and the search strategy are set.
                if gen < 0.2 * G
                    epsilon_k = 0;
                elseif gen < Tc
                    epsilon_k = (1 - tao) * epsilon_k;
                    if epsilon_k <= threshold
                        epsilon_k = (gen/G) * epsilon_0;
                    end
                end
                
                MatingPool1 = TournamentSelection(2,Problem.N,D_Dec,D_Pop,Fitness1);
                MatingPool2 = TournamentSelection(2,Problem.N,D,Fitness2);
                Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                [Population1,Fitness1,D_Dec,D_Pop] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N);
                [Population2,Fitness2,D] = EnvironmentalSelectionDec([Population2,Offspring1,Offspring2],Problem.N,epsilon_k);
                
                gen = ceil(Problem.FE/(2 * Problem.N));
                objs = sum(Population2.objs,2);
                max_obj = max(objs);
                epsilon_0 = min(max_obj,epsilon_0);
            end
        end
    end
end
classdef SGEA<ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <dynamic>
% Steady-state and generational evolutionary algorithm

%------------------------------- Reference --------------------------------
% S. Jiang and S. Yang, A steady-state and generational evolutionary
% algorithm for dynamic multiobjective optimization, IEEE Transactions on
% Evolutionary Computation, 2017, 21(1): 65-82.
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
            %% Generate random population
            Population = Problem.Initialization();
            [Population,Archive,Fitness] = EnvironmentalSelection(Population,length(Population));
            Centroid   = mean(Archive.decs,1);
            % Reset the number of saved populations (only for dynamic optimization)
            Algorithm.save = sign(Algorithm.save)*inf;
            % Archive for storing all populations before each change
            AllPop = [];	
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Changed(Problem,Population)
                    % Save the population before the change
                    AllPop = [AllPop,Population];
                    % React to the change
                    [Population,Archive,Centroid] = ChangeResponse(Problem,Population,Archive,Centroid);
                    Fitness = CalFitness(Population.objs);
                end
                Elitists = Population;
                for i = 1 : Problem.N
                    if rand < 0.5
                        MatingPool = TournamentSelection(2,2,Fitness);
                        Offspring  = OperatorGAhalf(Problem,Population(MatingPool));
                    else
                        MatingPool = TournamentSelection(2,1,Fitness);
                        SubPop     = Archive(randi(end));
                        Offspring  = OperatorGAhalf(Problem,[Population(MatingPool),SubPop]);
                    end
                    [Population,Archive,Fitness] = EnvironmentalSelection([Population,Offspring],Problem.N);
                end
                [Population,Archive,Fitness] = EnvironmentalSelection([Population,Elitists],Problem.N);
                if Problem.FE >= Problem.maxFE
                    % Return all populations
                    Population = [AllPop,Population];
                    [~,rank]   = sort(Population.adds(zeros(length(Population),1)));
                    Population = Population(rank);
                end
            end
        end
    end
end
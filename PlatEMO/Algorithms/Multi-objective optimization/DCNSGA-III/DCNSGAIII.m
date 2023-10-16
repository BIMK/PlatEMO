classdef DCNSGAIII < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained>
% Dynamic constrained NSGA-III
% cp --- 5 --- Decrease trend of the dynamic constraint boundary

%------------------------------- Reference --------------------------------
% R. Jiao, S. Zeng, C. Li, S. Yang, and Y. S. Ong, Handling constrained 
% many-objective optimization problems via problem transformation, IEEE 
% Transactions on Cybernetics, 2021, 51(10): 4834-4847.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            cp                    = Algorithm.ParameterSet(5);
            
            %% Generate the reference points and random population
            [Z,Problem.N]         = UniformPoint(Problem.N, Problem.M);
            Population            = Problem.Initialization();
            
            %% Calculate the initial dynamic constraint boundary
            [initialE, ~]         = max(max(0,Population.cons), [], 1);
            initialE(initialE==0) = 1;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Reduce the dynamic constraint boundry
                epsn       = ReduceBoundary(initialE, ceil(Problem.FE/Problem.N), ceil(Problem.maxFE/Problem.N)-1, cp);
                % Mating selection which prefers to select epsilon-feasible solutions
                MatingPool = TournamentSelection(2, Problem.N,sum(max(0, Population.cons-epsn), 2));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                % Update the reference points only consider the epsilon-feasible solution set
                Zmin       = min([Population(sum(max(0, Population.cons)<=epsn, 2)==size(Population.cons, 2)).objs; Offspring(sum(max(0,Offspring.cons)<=epsn, 2)==size(Offspring.cons, 2)).objs], [], 1);
                % Environment selection
                Population = EnvironmentalSelection([Population, Offspring], Problem.N, Z, Zmin, initialE, epsn);
            end
        end
    end
end
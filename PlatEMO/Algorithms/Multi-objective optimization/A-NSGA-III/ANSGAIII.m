classdef ANSGAIII < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Adaptive NSGA-III

%------------------------------- Reference --------------------------------
% H. Jain and K. Deb, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part II:
% Handling constraints and extending to an adaptive approach, IEEE
% Transactions on Evolutionary Computation, 2014, 18(4): 602-622.
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
            %% Generate the reference points and random population
            % All the reference points
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Z = sortrows(Z);
            % Distance between two consecutive reference points for the adaption
            interval = Z(1,end) - Z(2,end);
            % Initial population
            Population = Problem.Initialization();
            % Ideal point
            Zmin = min(Population(all(Population.cons<=0,2)).objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
                Z          = Adaptive(Population.objs,Z,Problem.N,interval);
            end
        end
    end
end
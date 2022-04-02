classdef AGEMOEA2 < ALGORITHM
    % <multi/many> <real/binary/permutation> <constrained/none>
    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
              MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
              Offspring  = OperatorGA(Population(MatingPool));
              [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end
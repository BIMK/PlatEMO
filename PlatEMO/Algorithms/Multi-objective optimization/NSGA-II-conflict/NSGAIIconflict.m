classdef NSGAIIconflict < ALGORITHM
% <many> <real/integer/label/binary/permutation>
% NSGA-II with conflict-based partitioning strategy
% NS     ---  2 --- Number of subspaces
% cycles --- 10 --- Number of cycles

%------------------------------- Reference --------------------------------
% A. L. Jaimes, C. A. Coello Coello, H. Aguirre, and K. Tanaka, Objective
% space partitioning using conflict information for solving many-objective
% problems, Information Sciences, 2014, 268: 305-327.
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
            [NS,cycles] = Algorithm.ParameterSet(2,10);
            Gc          = ceil(Problem.maxFE/Problem.N/cycles);

            %% Generate random population
            Psi        = {1:Problem.M};
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N,Psi);

            %% Optimization
            phase = true;
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N,Psi);
                if ~phase && mod(ceil(Problem.FE/Problem.N),Gc)/Gc < 0.3
                    % Change to the approximation phase
                    Psi   = {1:Problem.M};
                    phase = true;
                elseif phase && mod(ceil(Problem.FE/Problem.N),Gc)/Gc >= 0.3
                    % Change to the partitioning phase
                    Psi   = ConflictPartition(Population.objs,NS);
                    phase = false;
                end
            end
        end
    end
end
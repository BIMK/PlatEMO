classdef AENSGAII < ALGORITHM
% <2020> <multi> <real/integer/label/binary/permutation> <constrained/none> <dynamic>
% Autoencoding NSGA-II

%------------------------------- Reference --------------------------------
% L. Feng, W. Zhou, W. Liu, Y. S. Ong, and K. C. Tan. Solving dynamic
% multiobjective problem via autoencoding evolutionary search. IEEE
% Transations on Cybernetics, 2020, 52(5): 2649-2662.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            ChangeCount = 0;
            NDS = cell(1,100);
            % Reset the number of saved populations (only for dynamic optimization)
            Algorithm.save = sign(Algorithm.save)*inf;
            
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            % Archive for storing all populations before each change
            AllPop = [];
            NDS{1} = Population(FrontNo==1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Changed(Problem,Population)
                    ChangeCount        = ChangeCount+1;
                    NDS{ChangeCount+1} = Population(FrontNo==1);
                    AllPop             = [AllPop,Population];
                    % React to the change
                    curr_NDS = NDS{ChangeCount+1};
                    his_NDS  = NDS{ChangeCount};
                    [Population,FrontNo,CrowdDis] = AE_prediction(Problem,curr_NDS,his_NDS,Population,Problem.N);
                end
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
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
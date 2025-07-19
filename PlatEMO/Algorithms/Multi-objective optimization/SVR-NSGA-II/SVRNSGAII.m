classdef SVRNSGAII < ALGORITHM
% <2019> <multi> <real/integer/label/binary/permutation> <constrained/none> <dynamic>
% Support vector regression based NSGA-II

%------------------------------- Reference --------------------------------
% L. Cao, L. Xu, E. D. Goodman, C. Bao, and S. Zhu. Evolutionary dynamic
% multiobjective optimization assisted by a support vector regression
% predictor. IEEE Transactions on Evolutionary Computation, 2019, 24(2):
% 305-319.
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
            ChangeCount = Algorithm.ParameterSet(0);
            % Reset the number of saved populations (only for dynamic optimization)
            Algorithm.save = sign(Algorithm.save)*inf;
            
            %% Generate random population
            Population           = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            % Archive for storing all populations before each change
            AllPop = [];
            NDS    = cell(1,200);
            NDS{1} = Population;
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Changed(Problem,Population)
                    % Save the population before the change
                    ChangeCount        = ChangeCount+1;
                    NDS{ChangeCount+1} = Population;
                    AllPop = [AllPop,Population];
                    % React to the change
                    Population           = Reinitialization(Problem,NDS,ChangeCount);
                    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,length(Population));
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
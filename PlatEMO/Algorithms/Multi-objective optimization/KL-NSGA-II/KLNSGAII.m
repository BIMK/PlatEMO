classdef KLNSGAII < ALGORITHM
 % <2021> <multi> <real/integer/label/binary/permutation> <constrained/none> <dynamic>
 % Knowledge learning based NSGA-II

%------------------------------- Reference --------------------------------
% Q. Zhao, B. Yan, Y. Shi, and M. Middendorf. Evolutionary dynamic
% multiobjective optimization via learning from historical search process.
% IEEE Transactions on Cybernetics, 2021, 52(7): 6119-6130.
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
            % Initialize counters of SourceTemp, generation and change
            [j,tau,ChangeCount,zeta] = deal(1,0,0,0.2);
            % Reset the number of saved populations (only for dynamic optimization)
            Algorithm.save = sign(Algorithm.save)*inf;
            
            %% Generate random population
            Population           = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            % Archive for storing all populations before each change
            AllPop        = [];
            SourceTemp    = cell(1,10); % for save source data in each time step
            SourceTemp{1} = Population;
            FinalPopDec   = cell(1,100); % for save final solutions of each time step
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Changed(Problem,Population)
                    % Save the population before the change
                    ChangeCount   = ChangeCount+1;
                    SourceTemp{1} = Population;
                    j      = 0;
                    FinalPopDec{ChangeCount} = Population;
                    AllPop = [AllPop,Population];
                    % React to the change
                    [Population,FrontNo,CrowdDis] = Reinitialization(Problem,Population,zeta);
                end  
                    MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                    j = j+1; 
                if tau>50
                    Offspring      = OperatorGA(Problem,Population(MatingPool));
                    OffspringLearn = DA_TwoLayer(SourceTemp{end-1}.decs,SourceTemp{end}.decs,Population.decs,Problem,Problem.D);
                    Offspring      = [Offspring,OffspringLearn];
                else
                    Offspring = OperatorGA(Problem,Population(MatingPool));
                end 
                    [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
                    SourceTemp{j} = Population;
                    tau = tau+1;
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
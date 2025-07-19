classdef SSIORL < ALGORITHM
% <2025> <single> <real/integer> <large/none> <constrained/none>
% Search space independent operator based deep reinforcement learning
% fileName       --- 'Weight.mat' --- File name for saving/loading weights
% Mode           ---            1 --- 1. Test 2. Training
% TrainN         ---          100 --- Population size of GA for training
% TrainOffspring ---          800 --- Offspring size of GA for training
% TrainFE        ---        80000 --- Maximum evaluations of GA for training

%------------------------------- Reference --------------------------------
% Y. Tian, Y. Liu, S. Yang, and X. Zhang. Deep reinforcement learning based
% on search space independent operators for black-box continuous
% optimization. IEEE/CAA Journal of Automatica Sinica, 2025.
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
            [fileName,Mode,TrainN,TrainOffspring,TrainFE] = Algorithm.ParameterSet('Weight.mat',1,100,800,80000);
            if Mode == 1    % Test mode
                % Load weights
                if ischar(fileName)
                    try
                        load(fullfile(fileparts(mfilename('fullpath')),fileName),'-mat','agent');
                    catch
                        error('Fail to load agent from %s. This algorithm should be trained before used.',fileName);
                    end
                end
                Population = Problem.Initialization();
                while Algorithm.NotTerminated(Population)
                   State      = CalState(Population.decs, FitnessSingle(Population), Problem);
                   Action     = getAction(agent,State);
                   Action     = Action{:};
                   Parent1    = TournamentSelection(ceil(max(1,Problem.N*(Action(1)+1)/2)), Problem.N/2, FitnessSingle(Population));
                   Parent2    = TournamentSelection(1,Problem.N/2,FitnessSingle(Population));   
                   Parent3    = TournamentSelection(1,Problem.N/2,FitnessSingle(Population));
                   Offspring  = OperatorPPO(Population ,[Parent1, Parent2, Parent3], Action, Problem);
                   Population = [Population,Offspring];
                   [~,rank]   = sort(FitnessSingle(Population));
                   Population = Population(rank(1:Problem.N));
                end
            else            % Training mode
                mainPPO(TrainN, TrainFE, TrainOffspring, Problem, fileName);
            end
        end
    end
end
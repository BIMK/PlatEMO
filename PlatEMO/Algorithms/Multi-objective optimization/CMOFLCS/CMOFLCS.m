classdef CMOFLCS < ALGORITHM
% <2026> <multi> <real> <constrained>
% Constrained multi-objective optimization with fuzzy logic control system

%------------------------------- Reference --------------------------------
% Y. Wang, C. Hu, X. Chu, X. Wu, W. Gong, and X. Yan. Fuzzy logic control
% system assisted operator selection for constrained multi-objective
% optimization. IEEE Transactions on Fuzzy Systems, 2026, 34(6): 1986-2000.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Generate the sampling points and random population
            [W,Problem.N]   = UniformPoint(Problem.N,Problem.M);
            numOperators    = 2;
            Rewards         = 1e-3*ones(numOperators,1);
            rewardResetFE   = 20000;
            maxChange       = 1;
            changeThreshold = 0.01;
            historyLength   = 20;
            Population1     = Problem.Initialization();
            Population2     = Problem.Initialization();
            Fitness1        = CalFitness(Population1.objs,Population1.cons);
            Fitness2        = CalFitness(Population2.objs);
            stage           = 1;
            Z               = min([Population2.objs;Population1.objs],[],1);
            IdealHistory    = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            NadirHistory    = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            MeanHistory     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            IdealPoint      = min([Population2.objs],[],1);
            MeanPoint       = sum(Population2.objs,1)/(Problem.N*Problem.M);
            NadirPoint      = max([Population2.objs],[],1);
            problemType     = 0;
            Archive         = ArchiveUpdate([Population1 Population2],Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population2)
                generation                 = ceil(Problem.FE/Problem.N);
                IdealHistory(generation,:) = IdealPoint;
                NadirHistory(generation,:) = NadirPoint;
                MeanHistory(generation,:)  = MeanPoint;
                if generation >= historyLength
                    maxChange = MaximumChange(IdealHistory,NadirHistory,MeanHistory,generation,historyLength);
                end
                if maxChange <= changeThreshold && stage == 1
                    stage = -1;
                    fr    = length(find(sum(max(Population2.cons,0),2)<=0))/Problem.N;
                    if fr <= 0
                        problemType = 4;
                    end
                end
                if mod(Problem.FE, rewardResetFE) == 0
                    Rewards       = 1e-3*ones(numOperators,1);
                    Probabilities = 1 / numOperators * ones(1, numOperators);
                else
                    Convergence   = ConvergenceIndicator(Population1);
                    Diversity     = DiversityIndicator(Population1);
                    RewardGA      = Rewards(1)/sum(Rewards);
                    RewardDE      = Rewards(2)/sum(Rewards);
                    Probabilities = FLCS(RewardGA, RewardDE, Diversity, Convergence);
                    Probabilities = Probabilities/sum(Probabilities);
                end
                r        = rand;
                operator = numOperators;
                for t = 1 : numOperators
                    if r <= sum(Probabilities(1:t))
                        operator = t;
                        break;
                    end
                end
                if operator == 1
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                    Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                    Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                else
                    MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                    Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                    Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                end
                Z = min([Z;Offspring1.objs;Offspring2.objs],[],1);
                [Population1,Fitness1,Next] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,true);
                if stage == 1 || problemType ~= 4
                    [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Problem.N,false);
                else
                    [Population2,Fitness2] = AngleConstrainedSelection([Population2,Offspring1,Offspring2],Problem.N,Z,W);
                end
                Rewards(operator) = Rewards(operator) + sum(Next(Problem.N+1:end))/length([Offspring1 Offspring2]);
                Z                 = min([Population2.objs;Population1.objs;Z],[],1);
                IdealPoint        = min([Population2.objs;IdealPoint],[],1);
                MeanPoint         = sum([Population2.objs;MeanPoint],1)/(Problem.N*Problem.M);
                NadirPoint        = max([Population2.objs;NadirPoint],[],1);
                Archive           = ArchiveUpdate([Archive Population1 Population2],Problem.N);
                if Problem.FE >= Problem.maxFE
                    Population2 = Archive;
                end
            end
        end
    end
end

function Convergence = ConvergenceIndicator(Population)
    N          = size(Population, 2);
    distances  = zeros(N, 1);
    idealPoint = min([Population.objs], [], 1);
    for i = 1 : N
        distance     = norm(Population(i).obj - idealPoint);
        distances(i) = distance;
    end
    D_avg = mean(distances);
    D_min = min(distances);
    D_max = max(distances);
    Convergence = (D_max - D_avg) / (D_max - D_min);
end

function Diversity = DiversityIndicator(Population)
    N         = size(Population, 2);
    distances = zeros(1,N*(N-1)/2);
    pair      = 0;
    for i = 1 : N-1
        for j = i+1 : N
            distance = norm(Population(i).obj - Population(j).obj);
            pair     = pair + 1;
            distances(pair) = distance;
        end
    end
    D_avg = mean(distances);
    D_min = min(distances);
    D_max = max(distances);
    Diversity = (D_avg -  D_min ) / (D_max - D_min);
end

function maxChange = MaximumChange(IdealHistory,NadirHistory,MeanHistory,generation,historyLength)
    delta_value = 1e-6 * ones(1,size(IdealHistory,2));
    rz          = abs((IdealHistory(generation,:) - IdealHistory(generation - historyLength + 1,:)) ./ max(IdealHistory(generation - historyLength + 1,:),delta_value));
    nrz         = abs((NadirHistory(generation,:) - NadirHistory(generation - historyLength + 1,:)) ./ max(NadirHistory(generation - historyLength + 1,:),delta_value));
    arz         = abs((MeanHistory(generation,:) - MeanHistory(generation - historyLength + 1,:)) ./ max(MeanHistory(generation - historyLength + 1,:),delta_value));
    maxChange   = max([rz, nrz,arz]);
end
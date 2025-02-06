classdef MMEAPSL < ALGORITHM
% <2024> <multi> <real/integer/label/binary/permutation> <multimodal>
% Multimodal multi-objective evolutionary algorithm assisted by Pareto set learning

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, and Y. Jin. Growing neural gas network-based
% surrogate-assisted Pareto set learning for multimodal multi-objective
% optimization. Swarm and Evolutionary Computation, 2024, 87: 101541.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Parameters for GNG network
            params.N = Problem.N;
            params.MaxIt = 50;
            params.L = 30;
            params.epsilon_b = 0.2;
            params.epsilon_n = 0.006;
            params.alpha = 0.5;
            params.delta = 0.995;
            params.T = 30;
            genFlag = [];
            MaxGen = ceil(Problem.maxFE/Problem.N);
            netInitialized = 0;
            Pop1get = 0;
            gen = 0;
            
            %% Generate random population
            Population = Problem.Initialization();
            
            %% calculate fitness of populations
            [Fitness,D_Dec,~] = CalFitness(Population.objs,Population.decs);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen = gen + 1;
                %% Initial the GNG network when can
                if ~netInitialized
                    NDNum = sum(Fitness<1);
                    if NDNum >= 2
                        net = InitilizeGrowingGasNet(Population,params,Fitness);
                        netInitialized = 1;
                    end
                end
                
                if ~netInitialized || gen < 0.2 * MaxGen
                    MatingPool = TournamentSelection(2,Problem.N,D_Dec,Fitness);
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                    [Population,Fitness,D_Dec,~,net,genFlag] = EnvironmentalSelectionOrig([Population,Offspring],Problem.N,Problem,params,net,genFlag);
                else
                    if Pop1get == 0
                        Population1 = Population;
                        Fitness1 = CalFitnessSup(Population1.decs,net.w);
                        Pop1get = 1;
                    end
                     % use the network in generating offspring
                    MatingPool1 = TournamentSelection(2,Problem.N,D_Dec,Fitness);
                    Offspring1  = OperatorGAhalf(Problem,Population(MatingPool1));
                    V = net.w;
                    MatingPool2 = randi(length(V),1,Problem.N);
                    Offspring2 = Problem.Evaluation(OperatorGAhalf(Problem,V(MatingPool2,:)));
                    
                    MatingPool1 = TournamentSelection(2,Problem.N,-Fitness1);
                    Offspring3  = OperatorGAhalf(Problem,Population1(MatingPool1));
                    
                    Offspring = [Offspring1,Offspring2,Offspring3];
                    
                    [Population,Fitness,D_Dec,~,net,genFlag] = EnvironmentalSelectionOrig([Population,Offspring],Problem.N,Problem,params,net,genFlag);
                    [Population1,Fitness1] = EnvironmentalSelectionSup([Population1,Offspring],Problem.N,net);
                end
            end
        end
    end
end
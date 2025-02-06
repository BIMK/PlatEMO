classdef AGSEA < ALGORITHM
% <2024> <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% Automated guiding vector selection-based evolutionary algorithm

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, and X. Zhang. Deep reinforcement learning assisted
% automated guiding vector selection for large-scale sparse multi-objective
% optimization. Swarm and Evolutionary Computation, 2024, 88: 101606.
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
            %% Population initialization
            % Calculate the fitness of each decision variable        
            [Fitness1,TDec,TMask,TempPop]      = FitnessCal(Problem,5);
            [Population,Dec,Mask,FitnessSpea2] = EnvironmentalSelection(TempPop,TDec,TMask,Problem.N); 
            num_feature = 14;
            max_act     = 3;
            inputn      = zeros(num_feature,1);
            outputn     = zeros(1,1);
            net         = newff(inputn,outputn,[10 10 10],{'tansig','purelin'},'trainlm');
            Memory      = [];
            action      = 1;
            Fitness3    = zeros(1,Problem.D);
            Memory      = UpdateMemory(Problem,Memory,action,Population,Mask,Population,Mask);
            clear TempPop TDec TMask;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)                
                MatingPool     = TournamentSelection(2,2*Problem.N,FitnessSpea2);
                LastPopulation = Population;
                LastMask = Mask;
                [action, Fitness,Fitness3] = UsingNet(Problem,Fitness1,net,Memory,num_feature,max_act,Mask,action,Fitness3);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness,Mask,num_feature,Memory);
                Offspring        = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FitnessSpea2] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
                Memory = UpdateMemory(Problem,Memory,action,LastPopulation,LastMask,Population,Mask);
                net    = TrainNet(Problem,net,Memory,num_feature,max_act);                
            end
        end
    end
end
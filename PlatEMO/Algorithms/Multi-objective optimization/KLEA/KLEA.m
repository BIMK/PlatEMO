classdef KLEA < ALGORITHM
% <2025> <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% Knowledge learning-based evolutionary algorithm

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, Y. Zhang, and X. Zhang. Knowledge learning-based
% dimensionality reduction for solving large-scale sparse multiobjective
% optimization problems. IEEE Transactions on Cybernetics, 2025.
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
            [Fitness,TDec,TMask,TempPop]       = FitnessCal(Problem,5);
            [Population,Dec,Mask,FitnessSpea2] = EnvironmentalSelection(TempPop,TDec,TMask,Problem.N);          
            max_act = 3;
            inputn  = zeros(Problem.D + 1,1);
            outputn = zeros(1,1);
            net     = newff(inputn,outputn,[10 10 10],{'tansig','purelin'},'trainlm');
            Memory  = [];
            action  = 1;           
            Memory  = UpdateMemory(Memory,action,Population,Mask,Population,Mask);
            clear TempPop TDec TMask;

            %% Optimization
            while Algorithm.NotTerminated(Population)                
                MatingPool     = TournamentSelection(2,2*Problem.N,FitnessSpea2);
                LastPopulation = Population;
                LastMask       = Mask;                
                [OffDec,OffMask,action]            = Operator(action,Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness,Mask,Memory);
                Offspring                          = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FitnessSpea2] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
                Memory = UpdateMemory(Memory,action,LastPopulation,LastMask,Population,Mask);
                action = UsingNet(Problem,net,Memory,max_act);
                net    = TrainNet(Problem,net,Memory,max_act);
            end
        end
    end
end
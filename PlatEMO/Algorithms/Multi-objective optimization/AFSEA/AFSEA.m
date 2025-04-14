classdef AFSEA < ALGORITHM
% <2024> <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% Adjoint feature-selection-based evolutionary algorithm

%------------------------------- Reference --------------------------------
% P. Zhang, H. Yin, Y. Tian, and X. Zhang. An adjoint feature-selection-
% based evolutionary algorithm for sparse large-scale multiobjective
% optimization. Complex & Intelligent Systems, 2024.
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
            TDec    = [];
            TMask   = [];
            TempPop = [];
            Fitness = zeros(1,Problem.D+1);
            for i = 1 : 1+4*any(Problem.encoding~=4)  
                standard_Mask = zeros(1,Problem.D);
                standard_Dec  = rand(1,Problem.D);
                standard      = standard_Mask.*standard_Dec;
                Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                Dec(:,Problem.encoding==4) = 1;
                Mask        = eye(Problem.D);
                TPopulation = Problem.Evaluation([standard;Dec.*Mask]);
                Population  = TPopulation(2:end);
                TDec        = [TDec;Dec];
                TMask       = [TMask;Mask];
                TempPop     = [TempPop,Population];
                Fitness     = Fitness + NDSort([TPopulation.objs,TPopulation.cons],inf);       
            end
            % Pick the best subsequence based on NDS-FS
            Standard_0 = Fitness(1);
            [bestScore,bestindex] = min(Fitness(2:end));
            L = 1;
            SequenceSet = num2cell(setdiff(find(Fitness),1));
            while bestScore <= Standard_0*L
                Score = [];
                bestSequence = SequenceSet{bestindex};
                SequenceSet  = num2cell(setdiff(find(Fitness),[1,bestSequence]));
                for i = 1 : length(SequenceSet)
                    SequenceSet(i) = {[bestSequence,SequenceSet{i}]};
                    Score(i)       = sum(Fitness(SequenceSet{i}));
                    i = i + 1;
                end
                [bestScore,bestindex] = min(Score);
                L = L + 1;
            end
            bestSequence = bestSequence - 1;
            Mask         = false(Problem.N,Problem.D);
            for i = 1 : Problem.N
                if rand() < 0.5
                    Mask(i,bestSequence) = 1;
                else
                    Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness(2:end))) = 1;
                end
            end
            % Generate initial population
            Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            Dec(:,Problem.encoding==4) = 1;
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N);
            Delta = zeros(1,Problem.D);

           %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                delta      = Relief(Problem,Population,FrontNo);
                Delta      = Delta + abs(delta);
                if all(delta==0)
                    [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness);
                else
                    [OffMask,OffDec] = OperatorDelta(Problem,Mask(MatingPool,:),Dec(MatingPool,:),Delta);
                end
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
            end
        end
    end
end
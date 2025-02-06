classdef LRMOEA < ALGORITHM
% <2024> <multi> <real/integer/binary> <large/none> <constrained/none> <sparse> <robust>
% Large-scale robust multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, L. Zhang, K. C. Tan, and X. Zhang. An evolutionary
% algorithm for solving large-scale robust multi-objective optimization
% problems. IEEE Transactions on Evolutionary Computation, 2024.
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
            [thea,eta] = Algorithm.ParameterSet(0.2,1.15);
            W = UniformPoint(Problem.N,Problem.M);
            if Problem.encoding(1)==4
                Dec = ones(Problem.N,Problem.D);
            else
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            end
            Mask = false(size(Dec));
            for i = 1 : Problem.N
                Mask(i,randperm(end,ceil(rand.^2*end))) = true;
            end
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N);
            obj =Population.objs;
            score    = ones(1,Problem.D);
            score    = FitnessCal(Mask,FrontNo,score);
            Arch     = archives(obj,Dec,Mask);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool       = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                [OffDec,OffMask] = LROperator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),score);
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
                score = sum(Arch.masks)/size(Arch,2);
                if size(score,2) == 1
                    score = ones(1,Problem.D);
                end
                TemArch = Initialization(Population,Dec,Mask,FrontNo);
                Arch    = ArchUpdate(Problem,Arch,TemArch,thea,eta);
                if Problem.FE >= Problem.maxFE
                    Population = Final(Problem,Arch,Problem.N,W,score);
                end
            end
        end
    end
end
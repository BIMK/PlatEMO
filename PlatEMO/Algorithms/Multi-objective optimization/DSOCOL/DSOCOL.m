classdef DSOCOL < ALGORITHM
% <2026> <multi> <real> <large/none> <constrained>
% Dual-swarm optimizer with collaborative orthogonal learning

%------------------------------- Reference --------------------------------
% Y. Wang, C. Hu, X. Wu, T. Zhang, W. Gong, and X. Yan. Collaborative
% orthogonal learning for constrained multi-objective optimization. IEEE
% Transactions on Evolutionary Computation, 2026.
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
            [W,~]        = UniformPoint(5,Problem.M);
            generation   = 0;
            Swarm1       = Problem.Initialization();
            Archive      = ArchiveUpdate(Swarm1,Problem.N);
            Swarm2       = Problem.Initialization();
            cons         = [Swarm1.cons;Swarm2.cons];
            cons(cons<0) = 0;
            epsilon0     = max(sum(cons,2))*size(Swarm1.cons,2);
            if epsilon0 <= 50
                epsilon0 = 200;
            end
            progress = 0;
            Fitness1 = CalFitness(Swarm1.objs,Swarm1.cons,epsilon0);
            Fitness2 = CalFitness(Swarm1.objs,Swarm1.cons,0);

            %% Optimization
            while Algorithm.NotTerminated(Swarm1)
                cp         = (-log(epsilon0)-6)/log(1-0.5);
                epsilon    = epsilon0*(1-progress)^cp;
                generation = generation+1;
                if generation==1 || mod(generation,75)==0
                    [TrendOffspring,OrthogonalOffspring] = CollaborativeLearning(Problem,Swarm1,Swarm2,W,epsilon);
                    [Swarm1,Fitness1] = EnvironmentalSelection([Swarm1,TrendOffspring,OrthogonalOffspring],Problem.N,epsilon);
                    [Swarm2,Fitness2] = NicheSubsetSelection([Swarm2,TrendOffspring,OrthogonalOffspring],Problem.N,W,0);
                end
                if length(Swarm1) >= 2
                    Rank = randperm(length(Swarm1),floor(length(Swarm1)/2)*2);
                else
                    Rank = [1,1];
                end
                Loser          = Rank(1:end/2);
                Winner         = Rank(end/2+1:end);
                Change         = Fitness1(Loser) <= Fitness1(Winner);
                Temp           = Winner(Change);
                Winner(Change) = Loser(Change);
                Loser(Change)  = Temp;
                Offspring1     = CompetitiveUpdate(Problem,Swarm1(Loser),Swarm1(Winner),Fitness1(Winner),1);
                if length(Swarm2) >= 2
                    Rank = randperm(length(Swarm2),floor(length(Swarm2)/2)*2);
                else
                    Rank = [1,1];
                end
                Loser             = Rank(1:end/2);
                Winner            = Rank(end/2+1:end);
                Change            = Fitness2(Loser) <= Fitness2(Winner);
                Temp              = Winner(Change);
                Winner(Change)    = Loser(Change);
                Loser(Change)     = Temp;
                Offspring2        = CompetitiveUpdate(Problem,Swarm2(Loser),Swarm2(Winner),Fitness2(Winner),2);
                Offspring         = [Offspring1,Offspring2];
                [Swarm1,Fitness1] = EnvironmentalSelection([Swarm1,Offspring],Problem.N,epsilon);
                [Swarm2,Fitness2] = NicheSubsetSelection([Swarm2,Offspring],Problem.N,W,0);
                Archive           = ArchiveUpdate([Swarm1,Swarm2,Archive],Problem.N);
                progress          = progress+1/(Problem.maxFE/(Problem.N));
                if Problem.FE >= Problem.maxFE
                    Swarm1 = Archive;
                end
            end
        end
    end
end
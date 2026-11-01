classdef UICMO < ALGORITHM
% <2025> <multi> <real> <constrained>
% Utilizing feasible non-dominated solution information for constrained multi-objective optimization

%------------------------------- Reference --------------------------------
% Y. Wang, C. Hu, X. Wu, Z. Zhou, X. Yan, and W. Gong. Utilizing feasible
% non-dominated solution information for constrained multi-objective
% optimization. Information Sciences, 2025, 699: 121812.
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
            if Problem.M == 3 && Problem.N ~= 91
                Problem.N = 91;
            end
            Population1  = Problem.Initialization();
            Distance     = pdist2(Population1.objs,Population1.objs);
            beta         = max(max(Distance));
            Fitness1     = CalFitness(Population1.objs,Population1.cons,beta);
            delta        = 0.9;
            [W,~]        = UniformPoint(Problem.N,Problem.M);
            T            = ceil(Problem.N/10);
            B            = pdist2(W,W);
            [~,B]        = sort(B,2);
            B            = B(:,1:T);
            Population2  = Problem.Initialization();
            Z            = min([Population2.objs;Population1.objs],[],1);
            Archive      = ArchiveUpdate([Population1,Population2],Problem.N,0);
            ReferenceSet = FeasibleReferenceSet(Archive);

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                generation       = ceil(Problem.FE/(2*Problem.N));
                maxGeneration    = ceil(Problem.maxFE/(2*Problem.N));
                beta             = beta*(exp(-25*(generation/maxGeneration-0.5))/(1+exp(-25*(generation/maxGeneration-0.5))));
                MatingPool1      = TournamentSelection(2,Problem.N,Fitness1);
                Offspring1       = OperatorGA(Problem,Population1(MatingPool1));
                Z                = min([Z;Offspring1.objs],[],1);
                Offspring2       = [];
                Angle            = acos(1-pdist2(W,Offspring1.objs,'cosine'));
                [~,AssociateNew] = min(Angle,[],2);
                AssociatedOff    = Offspring1(AssociateNew);
                for i = 1 : Problem.N
                    if rand < delta
                        P = B(i,randperm(size(B,2)));
                    else
                        P = randperm(Problem.N);
                    end
                    Parent1     = FDPSelection(Population2(i),AssociatedOff(i),W,Z,ReferenceSet,i);
                    Parent2     = FDPSelection(Population2(P(1)),AssociatedOff(P(1)),W,Z,ReferenceSet,P(1));
                    Parent3     = FDPSelection(Population2(P(2)),AssociatedOff(P(2)),W,Z,ReferenceSet,P(2));
                    OffspringDE = OperatorDE(Problem,Parent1,Parent2,Parent3);

                    Offspring2     = [Offspring2,OffspringDE];
                    Population2(P) = FDPSelection(Population2(P),OffspringDE,W,Z,ReferenceSet,P);

                end
                Z            = min([Z;Offspring2.objs],[],1);
                Archive      = ArchiveUpdate([Archive,Offspring1,Offspring2],Problem.N,0);
                ReferenceSet = FeasibleReferenceSet(Archive);
                [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,beta);
                if Problem.FE >= Problem.maxFE
                    Population1 = Archive;
                end
            end
        end
    end
end

function Population = FeasibleReferenceSet(Population)
    fIndex     = all(Population.cons <= 0,2);
    Population = Population(fIndex);
    if isempty(Population)
        return
    else
        Fitness    = CalFitness(Population.objs,Population.cons,0);
        Next       = Fitness < 1;
        Population = Population(Next);

        N      = length(Population);
        PopObj = Population.objs;
        if N >= 3
            if size(PopObj,2) == 2
                [~,Extreme] = min(PopObj,[],1);
                [~,Rank]    = sort(PopObj(:,1));
                Reference   = [Extreme(1),Rank(ceil(N/2)),Extreme(2)];
            elseif size(PopObj,2) == 3
                [~,Reference] = min(PopObj,[],1);
            else
                return
            end
            Population = Population(unique(Reference,'stable'));
        end
    end
end
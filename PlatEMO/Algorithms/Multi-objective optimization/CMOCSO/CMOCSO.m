classdef CMOCSO < ALGORITHM
% <multi> <real> <large/none> <constrained>
% Competitive and cooperative swarm optimization constrained multi-objective optimization algorithm

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, D. Li, L. Wang, and L. Gao, A competitive and
% cooperative swarm optimizer for constrained multi-objective optimization
% problems, IEEE Transactions on Evolutionary Computation, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            CV         = sum(max(0,Population.cons),2);
            CVmax      = max(CV);
            epsilon_0  = CVmax;
            epsilon    = epsilon_0;
            Competitive_Pop = UpdateP1(Population,Problem.N,epsilon);
            [Cooperative_Pop,Cooperative_Pop_Fitness] = UpdateP2(Population,Problem.N);
            Population = UpdateP(Population,Problem.N);
            Tc    = 0.9 * ceil(Problem.maxFE/Problem.N);
            cp    = 2;
            alpha = 0.95;
            tao   = 0.05;
            y     = 10;
            G     = Problem.maxFE/Problem.N;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen       = ceil(Problem.FE/Problem.N);
                CV        = sum(max(0,Competitive_Pop.cons),2);
                CV_max    = max(CV);
                CVmax     = max([CV_max,CVmax]);
                epsilon_0 = CVmax;
                rf = sum(CV <= 1e-6) / length(Competitive_Pop);
                epsilon = update_epsilon(tao,epsilon,epsilon_0,rf,alpha,gen,Tc,cp);
                
                Competitive_Pop_Fitness = CalFitness(Competitive_Pop.objs,Competitive_Pop.cons,epsilon);
                
                if length(Competitive_Pop) >= 2
                    Rank = randperm(length(Competitive_Pop),floor(length(Competitive_Pop)/2)*2);
                else
                    Rank = [1,1];
                end
                Loser  = Rank(1:end/2);
                Winner = Rank(end/2+1:end);
                Change = Competitive_Pop_Fitness(Loser) <= Competitive_Pop_Fitness(Winner);
                Temp   = Winner(Change);
                Winner(Change) = Loser(Change);
                Loser(Change)  = Temp;
                Offspring1 = CompetitiveOperator(Problem,Competitive_Pop(Loser),Competitive_Pop(Winner),y);
                
                LearningPool = TournamentSelection(2,Problem.N,Cooperative_Pop_Fitness);
                Offspring2   = CooperativeOperator(Problem,Cooperative_Pop(LearningPool));
                
                Offspring = [Offspring1,Offspring2];
                
                Population = UpdateP([Population,Offspring],Problem.N);
                Competitive_Pop = UpdateP1([Competitive_Pop,Offspring],Problem.N,epsilon);
                gen = ceil(Problem.FE/Problem.N);
                y   = (Problem.M)^2*((gen/G)-1)^2+1;
                [Cooperative_Pop,Cooperative_Pop_Fitness] = UpdateP2([Offspring,Cooperative_Pop],Problem.N);
            end
        end
    end
end

function epsilon = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp)
% Update the value of epsilon

    if gen > Tc
        epsilon = 0;
    else
        if rf < alpha
            epsilon = (1 - tao) * epsilon_k;
        else
            epsilon = epsilon_0 * ((1 - (gen / Tc)) ^ cp);
        end
    end
end
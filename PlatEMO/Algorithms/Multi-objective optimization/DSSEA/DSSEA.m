classdef DSSEA < ALGORITHM
% <2025> <multi/many> <real/integer> <large> <constrained>
% Dynamic subspace search-based evolutionary algorithm

%------------------------------- Reference --------------------------------
% X. Ban, J. Liang, K. Yu, B. Qu, K. Qiao, P. N. Suganthan, and Y. Wang. A
% subspace search-based evolutionary algorithm for large-scale constrained
% multi-objective optimization and application. IEEE Transactions on
% Cybernetics, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [nSel,nPer] = Algorithm.ParameterSet(2,4);
    
            %% Generate random population
            Population = Problem.Initialization();
            Rank_DV   = ALDVA(Problem,Population,nSel,nPer);
    
            cons     = Population.cons;
            cons(cons<=0) = 0;
            conss    = sum(cons,2);
            epsilon0 = sum(conss);
    
            if epsilon0 == 0
                epsilon0 = 100;
            end
    
            Fitness = CalFitness_E(Population.objs,Population.cons,epsilon0);
            X = 0;
    
            %% Optimization
            while Algorithm.NotTerminated(Population)
                cp = (-log(epsilon0)-6)/log(1-0.75);
                epsilon = epsilon0*(1-X)^cp;
    
                IP = (Rank_DV-length(Rank_DV))/(1-length(Rank_DV));
                OP = (1-IP)*tanh(9*Problem.FE/Problem.maxFE)+IP;
    
                rand_num = rand(1,length(Rank_DV));
                DV_OP    = (rand_num<=OP);
                DV_NOP   = (rand_num>OP);
    
                if rand <= 0.5
                    Offspring = OperatorDE_pbest_1_main_OP(Population, Problem.N, Problem, Fitness, DV_OP, DV_NOP, 0.1);
                else
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness);
                    Offspring   = OperatorGA_OP(Problem,Population, Population(MatingPool1),DV_OP,DV_NOP);
                end
    
                [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Problem.N,epsilon,true);
    
                X = X + 1/(Problem.maxFE/Problem.N);
            end
        end
    end
end
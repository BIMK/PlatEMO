classdef  IMTCMO_BS < ALGORITHM
% <2024> <multi/many> <real/integer/label/binary/permutation> <large> <constrained/none>
% Improved evolutionary multitasking-based CMOEA with bidirectional sampling

%------------------------------- Reference --------------------------------
% K. Qiao, J. Liang, K. Yu, W. Guo, C. Yue, B. Qu, and P. N. Suganthan.
% Benchmark problems for large-scale constrained multi-objective
% optimization with baseline results. Swarm and Evolutionary Computation,
% 2024, 86: 101504.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao
% If you have any question, please email qiaokangjia@yeah.net

    methods
        function main(Algorithm,Problem)
            %% Parameter settings
            [Nw,Ns,g] = Algorithm.ParameterSet(10,30,50);
            
            %% Initialization
            Population1 = Problem.Initialization();
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Zmin1       = min(Population1.objs,[],1);
            
            Population2 = Problem.Initialization();
            Fitness2    = CalFitness(Population2.objs,Population2.cons);
            Zmin2       = min(Population2.objs,[],1);
            
            
            [V0,~]    = UniformPoint(Problem.N,Problem.M);
            [Vs0,L]   = UniformPoint(floor(Problem.N/5),Problem.M);
            [RefV,Vs] = deal(V0,Vs0);
            
            X    = 0;
            cons = [Population1.cons;Population2.cons];
            cons(cons<0) = 0;
            aaa  = sum(cons,2);
            VAR0 = max(aaa(~isinf(aaa)));
            if VAR0 == 0
                VAR0 = 1;
            end
            cnt = 0;
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                
                %% Udate the epsilon value
                cp=(-log(VAR0)-6)/log(1-0.5);
                VAR = VAR0*(1-X)^cp;
                cnt= cnt+1;
                if cnt==1 || mod(cnt,g)==0
                    temp_Population = Population1;
                    if rand > 0.5
                        [GuidingSolution,samplepop] = Diversity_DirectedSampling(Problem,temp_Population,Ns,Nw,RefV,VAR);
                    else
                        [GuidingSolution,samplepop] = Convergence_DirectedSampling(Problem,temp_Population,Ns,Nw,RefV,VAR);
                    end
                    [Population1,Fitness1] = Main_task_EnvironmentalSelection([Population1,samplepop],Problem.N,true);
                    [Population2,Fitness2] = Auxiliray_task_EnvironmentalSelection([Population2,samplepop], Problem.N,VAR);
                end
                %% Offspring generation
                MatingPool = [Population1(randsample(Problem.N,Problem.N))];
                [Mate1,Mate2,Mate3]  = Neighbor_Pairing_Strategy(MatingPool,Zmin1);
                Offspring1(1:Problem.N/2) = OperatorDE_rand_1(Problem,Mate1(1:Problem.N/2), Mate2(1:Problem.N/2), Mate3(1:Problem.N/2));
                Offspring1(1+Problem.N/2:Problem.N) = OperatorDE_pbest_1_main(Population1, Problem.N/2, Problem, Fitness1, 0.1);
                
                MatingPool = [Population2(randsample(Problem.N,Problem.N))];
                [Mate1,Mate2,Mate3]  = Neighbor_Pairing_Strategy(MatingPool,Zmin2);
                Offspring2(1:Problem.N/2) = OperatorDE_rand_1(Problem,Mate1(1:Problem.N/2),Mate2(1:Problem.N/2),Mate3(1:Problem.N/2));
                Offspring2(1+Problem.N/2:Problem.N) = OperatorDE_pbest_1_main(Population2, Problem.N/2, Problem, Fitness2, 0.1);
                
                Zmin1       = min([Zmin1;Offspring1.objs],[],1);
                Zmin2       = min([Zmin2;Offspring2.objs],[],1);
                
                [Population1,Fitness1] = Main_task_EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,true);
                [Population2,Fitness2] = Auxiliray_task_EnvironmentalSelection( [Population2,Offspring2,Offspring1], Problem.N,VAR);
                
                X=X+1/(Problem.maxFE/Problem.N);
            end
        end
    end
end
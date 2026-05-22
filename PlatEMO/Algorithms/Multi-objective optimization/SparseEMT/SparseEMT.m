classdef SparseEMT < ALGORITHM
% <2025> <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% Sparse evolutionary multitasking
% rmp    ---   1 --- random mating probability
% eval   --- 100 --- Number of evaluations for Knowledge transfer
% solver ---   1 --- 1.MSKEA 2.S-NSGA-II 3.MGCEA

%------------------------------- Reference --------------------------------
% J. Jiang, X. Fang, H. Wang, P. Tong, Z. Liu, B. Su, and F. Han. Turning
% sparse large-scale multiobjective optimization into evolutionary
% multitasking. IEEE Transactions on Evolutionary Computation, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jing Jiang and Xiang Fang (email: jingj0608@126.com)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [rmp,eval,solver] = Algorithm.ParameterSet(1,100,1);

            %% Initialization
            [Fitness,FitnessObj,Tasks,Dec,Mask] = InitialTasks(Problem,3,solver); 
            Population = SOLUTION_SparseEMT(Dec,Mask,Dec,Mask,Tasks,3,Problem);    
            
            %% Multitasking Optimization
            while Algorithm.NotTerminated(Population)
                switch solver
                    case 1
                        if ~exist('pv')
                            pv = Fitness;
                            sv = zeros(1,Problem.D);
                            Last_temp_num = 0;
                        end
                        Tasks = DividePop(Problem,Tasks,Population,solver,{sv,pv,Last_temp_num});
                        [Population,Tasks,sv,pv,Last_temp_num] = MSKEA_KT_Operator(Problem,Tasks,eval,rmp,sv,pv,Last_temp_num);
                    case 2
                        Tasks = DividePop(Problem,Tasks,Population,solver);
                        [Population,Tasks] = SNSGA2_KT_Operator(Problem,Tasks,eval);
                    case 3
                        REAL          = any(Problem.encoding==1);
                        if ~exist('FitnessOpt')
                            if REAL
                                FitnessOpt  = kmeans(sum(FitnessObj)',2);
                            else
                                FitnessOpt = kmeans(FitnessObj',2);
                            end
                            [SparseRate,FitnessOpt] = MGCEA_ClusterLabel(Fitness,FitnessOpt);
                            FitnessOpt  = FitnessOpt';
                            NearStage   = ceil(Problem.FE/(Problem.maxFE/10));
                            [FitnessLayer,LayerMax] = MGCEA_UpdateLayer(SparseRate,NearStage,FitnessOpt,Problem,[]);
                        end
                        Tasks = DividePop(Problem,Tasks,Population,solver,{Fitness,SparseRate,NearStage,FitnessLayer,LayerMax});
                        [Population,Tasks,FitnessOpt,NearStage,FitnessLayer,LayerMax] = MGCEA_KT_Operator(Problem,Tasks,eval,rmp,FitnessOpt,SparseRate,NearStage,FitnessLayer,LayerMax);
                end
            end
        end
    end
end
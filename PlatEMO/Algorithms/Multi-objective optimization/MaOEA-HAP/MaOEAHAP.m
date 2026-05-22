classdef MaOEAHAP < ALGORITHM
% <2026> <multi/many> <real/integer/label/binary/permutation>
% Hyper-curvature balanced indicator and adaptive phase exploration co-driven many-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% X. Yue, W. Wen, Y. Jiang, Y. Tian, and H. Peng. A hyper-curvature
% balanced indicator and adaptive phase exploration co-driven evolutionary
% algorithm for many-objective optimization. Swarm and Evolutionary
% Computation, 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            CAsize = Problem.N;
            p      = 1/Problem.M;
            k      = 0;

            %% Generate random population
            Population = Problem.Initialization();
            CA = UpdateCA([],Population,CAsize,p);
            DA = UpdateDA([],Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(DA)
                [ParentC,ParentM] = MatingSelection(CA,DA,Problem.N);
                Offspring         = [OperatorGA(Problem,ParentC,{1,20,0,0}),OperatorGA(Problem,ParentM,{0,0,1,20})];
                CA = UpdateCA(CA,Offspring,CAsize,p);
                DA = UpdateDA(DA,Offspring,Problem.N,p,Problem,k);
                k  = flog(DA,Problem.N,p);
            end
        end
    end
end
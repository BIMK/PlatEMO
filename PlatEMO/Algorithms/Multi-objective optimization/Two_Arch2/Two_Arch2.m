classdef Two_Arch2 < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Two-archive algorithm 2
% CAsize --- --- Convergence archive size
% p      --- --- The parameter of fractional distance

%------------------------------- Reference --------------------------------
% H. Wang, L. Jiao, and X. Yao, Two_Arch2: An improved two-archive
% algorithm for many-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2015, 19(4): 524-541.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [CAsize,p] = Algorithm.ParameterSet(Problem.N,1/Problem.M);

            %% Generate random population
            Population = Problem.Initialization();
            CA = UpdateCA([],Population,CAsize);
            DA = UpdateDA([],Population,Problem.N,p);

            %% Optimization
            while Algorithm.NotTerminated(DA)
                [ParentC,ParentM] = MatingSelection(CA,DA,Problem.N);
                Offspring         = [OperatorGA(Problem,ParentC,{1,20,0,0}),OperatorGA(Problem,ParentM,{0,0,1,20})];
                CA = UpdateCA(CA,Offspring,CAsize);
                DA = UpdateDA(DA,Offspring,Problem.N,p);
            end
        end
    end
end
classdef LMOEADS < ALGORITHM
% <multi> <real/integer> <large/none>
% Large-scale evolutionary multi-objective optimization assisted by
% directed sampling
% Nw --- 10 --- Cluster number
% Ns --- 30 --- The number of random sampling along each guiding direction

%------------------------------- Reference --------------------------------
% S. Qin, C. Sun, Y. Jin, Y. Tan, and J. Fieldsend, Large-scale
% evolutionary multi-objective optimization assisted by directed sampling,
% IEEE Transactions on Evolutionary Computation, 2021, 25(4): 724-738.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Shufen Qin
% E-mail: shufen.qin@stu.tyust.edu.cn

    methods
        function main(Algorithm,Problem)
            %% Parameter settings
            [Nw,Ns] = Algorithm.ParameterSet(10,30);

            %% Initialization
            Population = Problem.Initialization();
            [RefV,~]   = UniformPoint(Problem.N,Problem.M);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                GuidingSolution = DirectedSampling(Problem,Population,Ns,Nw,RefV);
                Population      = DoubleReproduction(Problem,Population,GuidingSolution,RefV);
            end
        end
    end
end
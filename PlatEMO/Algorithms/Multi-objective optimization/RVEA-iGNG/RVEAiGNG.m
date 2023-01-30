classdef RVEAiGNG < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% RVEA based on improved growing neural gas
% alpha --- 2 --- The parameter controlling the rate of change of penalty

%------------------------------- Reference --------------------------------
% Q. Liu, Y. Jin, M. Heiderich, T. Rodemann, and G. Yu, An adaptive
% reference vector-guided evolutionary algorithm using growing neural gas
% for many-objective optimization of irregular problems, IEEE Transactions
% on Cybernetics, 2022, 52(5): 2698-2711.
%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            params.N = Problem.N;
            params.MaxIt = 50;
            params.L = 50;
            params.epsilon_b = 0.2;
            params.epsilon_n = 0.006;
            params.alpha = 0.5;
            params.delta = 0.995;
            params.T = 50;
            
            %% Parameter setting
            alpha = Algorithm.ParameterSet(2);
            [V,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            net = InitilizeGrowingGasNet(V,Population,params);
            Archive = UpdateArchive(Population,[],Problem.N);
            scale = ones(1,Problem.M);
            zmin = min(Population.objs,[],1);
            genFlag = [];

            while Algorithm.NotTerminated(Population)
                MatingPool = randi(length(Population),1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                zmin       = min([zmin;Offspring.objs],[],1); 
                [Population,net,V,Archive,scale,genFlag] = EnvironmentalSelection([Population,Offspring],V,(Problem.FE/Problem.maxFE)^alpha,net,params,Archive,Problem,scale,zmin,genFlag);
            end
        end
    end
end
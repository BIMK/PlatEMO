function MOEADM2M(Global)
% <algorithm> <M>
% MOEA/D based on MOP to MOP
% K --- 10 --- Number of reference vectors

%------------------------------- Reference --------------------------------
% H. Liu, F. Gu, and Q. Zhang, Decomposition of a multiobjective
% optimization problem into a number of simple multiobjective subproblems,
% IEEE Transactions on Evolutionary Computation, 2014, 18(3): 450-455.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    K = Global.ParameterSet(10);

    %% Generate random population
    [W,K]      = UniformPoint(K,Global.M);
    Global.N   = ceil(Global.N/K)*K;
    S          = Global.N/K;
    Population = Global.Initialization();
    Population = Associate(Population,W,S);
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPoolLocal      = randi(S,S,K) + repmat(0:S:S*(K-1),S,1);
        MatingPoolGlobal     = randi(Global.N,1,Global.N);
        rnd                  = rand(S,K) < 0.7;
        MatingPoolLocal(rnd) = MatingPoolGlobal(rnd);
        Offspring  = Operator(Population,Population(MatingPoolLocal(:)));
        Population = Associate([Population,Offspring],W,S);
    end
end
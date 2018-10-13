function MOEADM2M(Global)
% <algorithm> <H-N>
% Decomposition of a Multiobjective Optimization Problem into a Number of
% Simple Multiobjective Subproblems
% K --- 10 --- Number of reference vectors
% operator --- MOEADM2M_operator

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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
        Offspring  = Global.Variation(Population([1:Global.N,MatingPoolLocal(:)']),Global.N,@MOEADM2M_operator);
        Population = Associate([Population,Offspring],W,S);
    end
end
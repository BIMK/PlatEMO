function tDEA(Global)
% <algorithm> <T>
% ¦È-dominance based evolutionary algorithm

%------------------------------- Reference --------------------------------
% Y. Yuan, H. Xu, B. Wang, and X. Yao, A new dominance relation-based
% evolutionary algorithm for many-objective optimization, IEEE Transactions
% on Evolutionary Computation, 2016, 20(1): 16-37.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference points and random population
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    Population   = Global.Initialization();
    [z,znad]     = deal(min(Population.objs),max(Population.objs));

    %% Optimization
    while Global.NotTermination(Population) 
        MatingPool = randi(Global.N,1,Global.N);
        Offspring  = GA(Population(MatingPool));
        [Population,z,znad] = EnvironmentalSelection([Population,Offspring],W,Global.N,z,znad);
    end
end
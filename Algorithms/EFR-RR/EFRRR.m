function EFRRR(Global)
% <algorithm> <E>
% Ensemble fitness ranking with a ranking restriction scheme
% K --- 2 --- Number of nearest weight vectors

%------------------------------- Reference --------------------------------
% Y. Yuan, H. Xu, B. Wang, B. Zhang, and X. Yao, Balancing convergence and
% diversity in decomposition-based many-objective optimizers, IEEE
% Transactions on Evolutionary Computation, 2016, 20(2): 180-198.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    K = Global.ParameterSet(2);

    %% Generate the reference points and random population
    [W,Global.N]    = UniformPoint(Global.N,Global.M);
    Population      = Global.Initialization();
    [PopObj,z,znad] = Normalization(Population.objs,min(Population.objs),max(Population.objs));
    RgFrontNO       = MaximumRanking(PopObj,W,K);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,RgFrontNO);
        Offspring  = GA(Population(MatingPool));
        [Population,z,znad] = EnvironmentalSelection([Population,Offspring],W,Global.N,K,z,znad);
    end
end
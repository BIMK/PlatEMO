function EFRRR(Global)
% <algorithm> <A-G>
% Balancing Convergence and Diversity in Decomposition-Based Many-Objective
% Optimizers
% K --- 2 --- Number of nearest weight vectors

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,z,znad] = EnvironmentalSelection([Population,Offspring],W,Global.N,K,z,znad);
    end
end
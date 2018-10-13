function RPDNSGAII(Global)
% <algorithm> <O-Z>
% A New Decomposition-Based NSGA-II for Many-Objective Optimization

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference points and random population
    [RPSet,Global.N] = UniformPoint(Global.N,Global.M);
    Population       = Global.Initialization();
    [~,FrontNo,d2]   = EnvironmentalSelection(Population,RPSet,Global.N);

    %% Optimization
    while Global.NotTermination(Population) 
        MatingPool = TournamentSelection(2,Global.N,FrontNo,d2);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,d2] = EnvironmentalSelection([Population,Offspring],RPSet,Global.N);
    end
end
function NSGAIISDR(Global)
% <algorithm> <H-N>
% A Strengthened Dominance Relation Considering Convergence and Diversity
% for Evolutionary Many-Objective Optimization

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
    zmin       = min(Population.objs,[],1);
    zmax       = max(Population.objs,[],1);
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N,zmin,zmax);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        zmin       = min([zmin;Offspring.objs],[],1);
        zmax       = max(Population(FrontNo==1).objs,[],1);
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,zmin,zmax);
    end
end
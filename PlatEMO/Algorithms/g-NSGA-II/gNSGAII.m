function gNSGAII(Global)
% <algorithm> <G>
% g-dominance based NSGA-II
% Point --- --- Preferred point

%------------------------------- Reference --------------------------------
% J. Molina, L. V. Santana, A .G. Hernandez-Diaz, C. A. Coello Coello, and
% R. Caballero, g-dominance: Reference point based dominance for
% multiobjective metaheuristics, European Journal of Operational Research,
% 2009, 197(2): 685-692.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Point = Global.ParameterSet(zeros(1,Global.M)+0.5);

    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NDSort(Evaluate(Population.objs,Point),inf);
    CrowdDis   = CrowdingDistance(Population.objs,FrontNo);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = GA(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,Point);
    end
end
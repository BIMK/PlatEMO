function CPSMOEA(Global)
% <algorithm> <A-G>
% A Classification and Pareto Domination based Multiobjective Evolutionary
% Algorithm
% M ---  3 --- Number of generated offsprings for each solution
% operator --- DE

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    M = Global.ParameterSet(3);

    %% Generate random population
    Population   = Global.Initialization();
    [Pgood,Pbad] = NDS(Population,floor(Global.N/2));

    %% Optimization
    while Global.NotTermination(Population)
        KNN(Pgood.decs,Pbad.decs);
        Offspring  = GenerateOffsprings(Global,Population,M);
        Population = NDS([Population,Offspring],Global.N);
        FrontNo    = NDSort(Offspring.objs,1);
        Pgood      = NDS([Pgood,Offspring(FrontNo==1)],floor(Global.N/2));
        Pbad       = NDS([Pbad,Offspring(FrontNo~=1)],floor(Global.N/2));
    end
end
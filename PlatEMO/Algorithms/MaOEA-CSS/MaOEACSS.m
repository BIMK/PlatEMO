function MaOEACSS(Global)
% <algorithm> <M>
% Many-objective evolutionary algorithms based on coordinated selection
% strategy
% t --- 0 --- Threshold value in environmental selection

%------------------------------- Reference --------------------------------
% Z. He and G. G. Yen, Many-objective evolutionary algorithms based on
% coordinated selection strategy, IEEE Transactions on Evolutionary
% Computation, 2017, 21(2): 220-233.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    t = Global.ParameterSet(0);

    %% Generate random population
    Population = Global.Initialization();
    Zmin       = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,Zmin);
        Offspring  = GA(Population(MatingPool));
        Zmin       = min([Zmin;Offspring.objs],[],1);
        Population = EnvironmentalSelection([Population,Offspring],Zmin,t,Global.N);
    end
end
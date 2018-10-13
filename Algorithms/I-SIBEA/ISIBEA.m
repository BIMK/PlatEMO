function ISIBEA(Global)
% <algorithm> <H-N>
% An interactive simple indicator-based evolutionary algorithm (I-SIBEA)
% for multiobjective optimization problems
% Point --- --- Preferred point

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Point = Global.ParameterSet(ones(1,Global.M));

    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NDSort(Evaluate(Population.objs,Point),inf);
    WHVLoss    = CalWHVLoss(Population.objs,FrontNo);
    wz = [];    % Weight distribution function value
    AA = [];    % Preferred set
    RA = [];    % Non-preferred set

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-WHVLoss);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,WHVLoss] = EnvironmentalSelection([Population,Offspring],Global.N,wz,AA,RA);
        if ~mod(Global.gen,ceil(Global.maxgen/4))
            [wz,AA,RA] = Interaction(Population.objs,Point);
        end
    end
end
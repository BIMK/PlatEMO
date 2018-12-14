function ISIBEA(Global)
% <algorithm> <I>
% Interactive simple indicator-based evolutionary algorithm
% Point --- --- Preferred point

%------------------------------- Reference --------------------------------
% T. Chugh, K. Sindhya, J. Hakanen, and K. Miettinen, An interactive simple
% indicator-based evolutionary algorithm (I-SIBEA) for multiobjective
% optimization problems, Proceedings of the International Conference on
% Evolutionary Multi-Criterion Optimization, 2015, 277-291.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
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
        Offspring  = GA(Population(MatingPool));
        [Population,FrontNo,WHVLoss] = EnvironmentalSelection([Population,Offspring],Global.N,wz,AA,RA);
        if ~mod(Global.gen,ceil(Global.maxgen/4))
            [wz,AA,RA] = Interaction(Population.objs,Point);
        end
    end
end
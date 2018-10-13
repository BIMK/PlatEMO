function KnEA(Global)
% <algorithm> <H-N>
% A Knee Point Driven Evolutionary Algorithm for Many-Objective
% Optimization
% rate --- 0.5 --- Rate of knee points in the population

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    rate = Global.ParameterSet(0.5);

    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NDSort(Population.objs,Population.cons,inf);
    KneePoints = zeros(1,Global.N);     % Set of knee points
    r          = -ones(1,2*Global.N);	% Ratio of size of neighorhood
    t          = -ones(1,2*Global.N);	% Ratio of knee points

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,FrontNo,KneePoints);
        Offspring  = Global.Variation(Population(MatingPool));
        Population = [Population,Offspring];
        [FrontNo,MaxFNo]                = NDSort(Population.objs,Population.cons,Global.N);
        [KneePoints,Distance,r,t]       = FindKneePoints(Population.objs,FrontNo,MaxFNo,r,t,rate);
        [Population,FrontNo,KneePoints] = EnvironmentalSelection(Population,FrontNo,MaxFNo,KneePoints,Distance,Global.N);      
    end
end
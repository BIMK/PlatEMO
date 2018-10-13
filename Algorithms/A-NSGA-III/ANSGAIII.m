function ANSGAIII(Global)
% <algorithm> <A-G>
% An Evolutionary Many-Objective Optimization Algorithm Using
% Reference-Point Based Nondominated Sorting Approach, Part II: Handling
% Constraints and Extending to an Adaptive Approach

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference points and random population
    % All the reference points
    [Z,Global.N] = UniformPoint(Global.N,Global.M);
    Z = sortrows(Z);
    % Distance between two consecutive reference points for the adaption
    interval = Z(1,end) - Z(2,end);
    % Initial population
    Population = Global.Initialization();
    % Ideal point
    Zmin = min(Population(all(Population.cons<=0,2)).objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,sum(max(0,Population.cons),2));
        Offspring  = Global.Variation(Population(MatingPool));
        Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
        Population = EnvironmentalSelection([Population,Offspring],Global.N,Z,Zmin);
        Z          = Adaptive(Population.objs,Z,Global.N,interval);
    end
end
function ARMOEA(Global)
% <algorithm> <A>
% Adaptive reference points based multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% Y. Tian, R. Cheng, X. Zhang, and Y. Jin, An indicator-based
% multiobjective evolutionary algorithm with reference point adaptation for
% better versatility, IEEE Transactions on Evolutionary Computation, 2018,
% 22(4): 609-622.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the sampling points and random population
    Population = Global.Initialization();
    W          = UniformPoint(Global.N,Global.M);
    [Archive,RefPoint,Range] = UpdateRefPoint(Population(all(Population.cons<=0,2)).objs,W,[]);
    
    %% Start the interations
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population,RefPoint,Range);
        Offspring  = GA(Population(MatingPool));
        [Archive,RefPoint,Range] = UpdateRefPoint([Archive;Offspring(all(Offspring.cons<=0,2)).objs],W,Range);
        [Population,Range]       = EnvironmentalSelection([Population,Offspring],RefPoint,Range,Global.N);
    end
end
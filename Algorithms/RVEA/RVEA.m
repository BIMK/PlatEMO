function RVEA(Global)
% <algorithm> <O-Z>
% A Reference Vector Guided Evolutionary Algorithm for Many-objective
% Optimization
% alpha ---   2 --- The parameter controlling the rate of change of penalty
% fr    --- 0.1 --- The frequency of employing reference vector adaptation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [alpha,fr] = Global.ParameterSet(2,0.1);

    %% Generate the reference points and random population
    [V0,Global.N] = UniformPoint(Global.N,Global.M);
    Population    = Global.Initialization();
    V             = V0;

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = randi(length(Population),1,Global.N);
        Offspring  = Global.Variation(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],V,(Global.gen/Global.maxgen)^alpha);
        if ~mod(Global.gen,ceil(fr*Global.maxgen))
            V(1:Global.N,:) = ReferenceVectorAdaptation(Population.objs,V0);
        end
    end
end
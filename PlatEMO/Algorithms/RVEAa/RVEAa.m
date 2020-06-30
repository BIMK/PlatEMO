function RVEAa(Global)
% <algorithm> <R>
% RVEA embedded with the reference vector regeneration strategy
% alpha ---   2 --- The parameter controlling the rate of change of penalty
% fr    --- 0.1 --- The frequency of employing reference vector adaptation

%------------------------------- Reference --------------------------------
% R. Cheng, Y. Jin, M. Olhofer, and B. Sendhoff, A reference vector guided
% evolutionary algorithm for many-objective optimization, IEEE Transactions
% on Evolutionary Computation, 2016, 20(5): 773-791.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [alpha,fr] = Global.ParameterSet(2,0.1);

    %% Generate the reference points and random population
    [V0,Global.N] = UniformPoint(Global.N,Global.M);
    Population    = Global.Initialization();
    V             = [V0;rand(Global.N,Global.M)];

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = randi(length(Population),1,Global.N);
        Offspring  = GA(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],V,(Global.gen/Global.maxgen)^alpha);
        if ~mod(Global.gen,ceil(fr*Global.maxgen))
            V(1:Global.N,:) = ReferenceVectorAdaptation(Population.objs,V0);
        end
        V(Global.N+1:end,:) = ReferenceVectorRegeneration(Population.objs,V(Global.N+1:end,:));
        if Global.evaluated >= Global.evaluation
            Population = Truncation(Population,Global.N);
        end
    end
end
function MyODEMR(Global)
% <algorithm> <M>
% Many-objective differential evolution with mutation restriction
% nP --- 500 --- Number of reference points for IGD calculation

%------------------------------- Reference --------------------------------
% R. Denysiuk, L. Costa, and I. E. Santo, Many-objective optimization using
% differential evolution with variable-wise mutation restriction,
% Proceedings of the 15th Annual Conference on Genetic and Evolutionary
% Computation, 2013, 591-598.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

    %% Parameter setting
    nP = Global.ParameterSet(500);

    %% Generate hyperplane 
    P = UniformPoint(nP,Global.M);

    %% Generate random population
    Population = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        Offspring  = Operator(Population(1:Global.N),Population(randi(Global.N,1,Global.N)),Population(randi(Global.N,1,Global.N)));
        Population = EnvironmentalSelection([Population,Offspring],Global.N,P);
    end
end
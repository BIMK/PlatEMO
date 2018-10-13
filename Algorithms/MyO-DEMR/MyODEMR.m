function MyODEMR(Global)
% <algorithm> <H-N>
% Many-Objective Optimization using Differential Evolution with
% Variable-Wise Mutation Restriction
% nP --- 500 --- Number of reference points for IGD calculation
% operator   --- DEMR

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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
        % Perform variation
        Offspring = Global.Variation(Population([1:Global.N,randi(Global.N,1,2*Global.N)]),Global.N,@DEMR);
        % Form new population
        [Population] = EnvironmentalSelection([Population,Offspring],Global.N,P);
    end
end
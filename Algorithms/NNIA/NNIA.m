function NNIA(Global)
% <algorithm> <H-N>
% Multiobjective Immune Algorithm with Nondominated Neighbor-Based
% Selection
% nA ---  20 --- Size of active population
% nC --- 100 --- Size of clone population

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [nA,nC] = Global.ParameterSet(20,100);

    %% Generate random population
    B = Global.Initialization();                % Antibody population
    D = UpdateDominantPopulation(B,Global.N);	% Dominant population

    %% Optimization
    while Global.NotTermination(D)
        A  = D(1:min(nA,length(D)));            % Active population
        C  = Cloning(A,nC);                     % Clone population
        C1 = Global.Variation([C,A(randi(length(A),1,length(C)))],length(C));
        D  = UpdateDominantPopulation([D,C1],Global.N);
    end
end
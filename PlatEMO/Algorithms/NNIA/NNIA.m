function NNIA(Global)
% <algorithm> <N>
% Nondominated neighbor immune algorithm
% nA ---  20 --- Size of active population
% nC --- 100 --- Size of clone population

%------------------------------- Reference --------------------------------
% M. Gong, L. Jiao, H. Du, and L. Bo, Multiobjective immune algorithm with
% nondominated neighbor-based selection, Evolutionary Computation, 2008,
% 16(2): 225-255.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
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
        C1 = GAhalf([C,A(randi(length(A),1,length(C)))]);
        D  = UpdateDominantPopulation([D,C1],Global.N);
    end
end
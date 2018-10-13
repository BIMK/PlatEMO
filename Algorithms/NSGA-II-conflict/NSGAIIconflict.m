function NSGAIIconflict(Global)
% <algorithm> <H-N>
% Objective space partitioning using conflict information for solving
% many-objective problems
% NS     ---  2 --- Number of subspaces
% cycles --- 10 --- Number of cycles

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [NS,cycles] = Global.ParameterSet(2,10);
    Gc          = ceil(Global.maxgen/cycles);

    %% Generate random population
    Psi        = {1:Global.M};
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N,Psi);

    %% Optimization
    phase = true;
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,Psi);
        if ~phase && mod(Global.gen,Gc)/Gc < 0.3
            % Change to the approximation phase
            Psi   = {1:Global.M};
            phase = true;
        elseif phase && mod(Global.gen,Gc)/Gc >= 0.3
            % Change to the partitioning phase
            Psi   = ConflictPartition(Population.objs,NS);
            phase = false;
        end
    end
end
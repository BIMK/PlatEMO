function BCEIBEA(Global)
% <algorithm> <A-G>
% Pareto or Non-Pareto: Bi-Criterion Evolution in Multi-Objective
% Optimization
% kappa --- 0.05 --- Fitness scaling factor

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    kappa = Global.ParameterSet(0.05);
    
    %% Generate random population
    NPC      = Global.Initialization();
    [PC,nND] = PCSelection(NPC,Global.N);
    
    %% Optimization
    while Global.NotTermination(PC)
        % PC evolving
        NewPC = Exploration(Global,PC,NPC,nND);
        
        % NPC selection
        NPC = EnvironmentalSelection([NPC,NewPC],Global.N,kappa);
        
        % NPC evolving
        MatingPool = TournamentSelection(2,Global.N,-CalFitness(NPC.objs,kappa));
        NewNPC     = Global.Variation(NPC(MatingPool));
        NPC        = EnvironmentalSelection([NPC,NewNPC],Global.N,kappa);
        
        % PC selection
        [PC,nND] = PCSelection([PC,NewNPC,NewPC],Global.N);
    end
end
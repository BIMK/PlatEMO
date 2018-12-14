function BCEIBEA(Global)
% <algorithm> <B>
% Bi-criterion evolution based IBEA
% kappa --- 0.05 --- Fitness scaling factor

%------------------------------- Reference --------------------------------
% M. Li, S. Yang, and X. Liu, Pareto or non-Pareto: Bi-criterion evolution
% in multiobjective optimization, IEEE Transactions on Evolutionary
% Computation, 2016, 20(5): 645-665.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
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
        NewPC = Exploration(PC,NPC,nND,Global.N);
        
        % NPC selection
        NPC = EnvironmentalSelection([NPC,NewPC],Global.N,kappa);
        
        % NPC evolving
        MatingPool = TournamentSelection(2,Global.N,-CalFitness(NPC.objs,kappa));
        NewNPC     = GA(NPC(MatingPool));
        NPC        = EnvironmentalSelection([NPC,NewNPC],Global.N,kappa);
        
        % PC selection
        [PC,nND] = PCSelection([PC,NewNPC,NewPC],Global.N);
    end
end
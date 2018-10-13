function MSOPSII(Global)
% <algorithm> <H-N>
% MSOPS-II: A general-purpose Many-Objective optimiser

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
    Feasible   = all(Population.cons<=0,2);
    if any(Feasible)
        Archive = UpdateArchive([],Population(Feasible),Global.N);
        Weight  = UpdateWeight([],Population(Feasible).objs,Global.N);
    else
        [~,best] = min(sum(max(0,Population.cons),2));
        Archive  = Population(best);
        Weight   = Population(best).objs;
    end

    %% Optimization
    % As the number of solutions in the archive is too large and
    % uncontrollable, use the population as the final output
    while Global.NotTermination(Population)
        Parents    = MatingSelection(Population,Archive);
        Offspring  = Global.Variation(Parents,Global.N);
        Feasible   = all(Offspring.cons<=0,2);
        Archive    = UpdateArchive(Archive,Offspring(Feasible),Global.N);
        Weight     = UpdateWeight(Weight,Offspring(Feasible).objs,Global.N);
        Population = EnvironmentalSelection([Population,Offspring],Weight,Global.N);
    end
end
function MSOPSII(Global)
% <algorithm> <M>
% Multiple single objective Pareto sampling II

%------------------------------- Reference --------------------------------
% E. J. Hughes, MSOPS-II: A general-purpose many-objective optimiser,
% Proceedings of the IEEE Congress on Evolutionary Computation, 2007,
% 3944-3951.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
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
        Offspring  = GAhalf(Parents);
        Feasible   = all(Offspring.cons<=0,2);
        Archive    = UpdateArchive(Archive,Offspring(Feasible),Global.N);
        Weight     = UpdateWeight(Weight,Offspring(Feasible).objs,Global.N);
        Population = EnvironmentalSelection([Population,Offspring],Weight,Global.N);
    end
end
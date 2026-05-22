function [state1, action1, nextstate1, ter, reward1] = GenerateSample(Problem, action,LastPopulation,Population)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    LastPopulationHV = HV(LastPopulation,max([LastPopulation.objs;Population.objs]));
    PopulationHV     = HV(Population,max([LastPopulation.objs;Population.objs]));
    reward           = (PopulationHV - LastPopulationHV)/LastPopulationHV;
    if isnan(reward)
        reward = 0;
    end

    LastObjsVar = var(LastPopulation.objs);
    LastObjsCon = sum(LastPopulation.objs);
    LastCV      = sum(max(Population.cons,0),'all');
    LastRatio   = Problem.FE/Problem.maxFE;
    LastState   = [LastObjsVar, LastObjsCon, LastCV, LastRatio];

    CurrentObjsVar = var(Population.objs);
    CurrentObjsCon = sum(Population.objs);
    CurrentCV      = sum(max(Population.cons,0),'all');
    CurrentRatio   = Problem.FE/Problem.maxFE;
    CurrentState   = [CurrentObjsVar, CurrentObjsCon, CurrentCV, CurrentRatio];

    state1     = LastState;
    action1    = action;
    nextstate1 = CurrentState;
    ter        = 0;
    reward1    = reward;
end
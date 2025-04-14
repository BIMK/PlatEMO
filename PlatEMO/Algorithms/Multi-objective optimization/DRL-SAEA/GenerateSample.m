function [state, action, nextstate, ter, reward] = GenerateSample(Problem, action,LastPopulation,Population)
    LastPopulationHV = HV(LastPopulation,max([LastPopulation.objs;Population.objs]));
    PopulationHV = HV(Population,max([LastPopulation.objs;Population.objs]));
    reward1 = (PopulationHV - LastPopulationHV)/LastPopulationHV;
    if isnan(reward1)
        reward1 = 0;
    end
    LastObjsVar = sum(var(LastPopulation.objs));
    LastObjsCon = sum(LastPopulation.objs,'all');
    LastCV = sum(max(0,LastPopulation.cons),'all');
    LastRatio = Problem.FE/Problem.maxFE;
    LastState = [LastObjsVar, LastObjsCon, LastCV, LastRatio];
    CurrentObjsVar = sum(var(Population.objs));
    CurrentObjsCon = sum(Population.objs,'all');
    CurrentCV = sum(max(0,Population.cons),'all');
    CurrentRatio = Problem.FE/Problem.maxFE;
    CurrentState = [CurrentObjsVar, CurrentObjsCon, CurrentCV, CurrentRatio];
    state = LastState;
    nextstate = CurrentState;
    ter = 0;
    if LastCV == 0
        reward2 = 0;        
    elseif CurrentCV < LastCV
        reward2 = abs((CurrentCV-LastCV)/LastCV);
    else
        reward2 = -abs((CurrentCV-LastCV)/LastCV);
    end  
    reward = reward1 + reward2;
    if action == -1
        action = 3;
    end
end
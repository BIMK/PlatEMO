function [state, action, nextstate, ter, reward] = GenerateSample(Problem, action,LastPopulation,Population)

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
    reward1          = (PopulationHV - LastPopulationHV)/LastPopulationHV;
    if isnan(reward1)
        reward1 = 0;
    end   
    LastPopDec     = LastPopulation.decs;
    LastSparse     = sum(LastPopDec~=0, 2) ./ size(LastPopDec,2);
    LastS1         = mean(LastSparse);
    LastS2         = var(LastSparse);
    LastObjsVar    = sum(var(LastPopulation.objs));
    LastObjsCon    = sum(LastPopulation.objs,'all');
    LastCV         = sum(max(0,LastPopulation.cons),'all');
    LastRatio      = sum(sum(max(0,LastPopulation.cons),2)==0)/Problem.N;
    LastState      = [LastS1, LastS2, LastObjsVar, LastObjsCon, LastCV, LastRatio];
    CurrentPopDec  = LastPopulation.decs;
    CurrentSparse  = sum(CurrentPopDec~=0, 2) ./ size(CurrentPopDec,2);
    CurrentS1      = mean(CurrentSparse);
    CurrentS2      = var(CurrentSparse);
    CurrentObjsVar = sum(var(Population.objs));
    CurrentObjsCon = sum(Population.objs,'all');
    CurrentCV      = sum(max(0,Population.cons),'all');
    CurrentRatio   = sum(sum(max(0,Population.cons),2)==0)/Problem.N;
    CurrentState   = [CurrentS1, CurrentS2,CurrentObjsVar, CurrentObjsCon, CurrentCV, CurrentRatio];
    state          = LastState;
    nextstate      = CurrentState;
    ter            = 0;
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
function Memory = UpdateMemory(Memory,action,LastPopulation,LastMask,Population,Mask)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    LastPopulationHV = HV(LastPopulation,max([LastPopulation.objs;Population.objs]));
    PopulationHV     = HV(Population,max([LastPopulation.objs;Population.objs]));
    reward           = (PopulationHV - LastPopulationHV)/LastPopulationHV;

    LastState    = sum(LastMask, 1)./size(LastMask, 1);
    CurrentState = sum(Mask, 1)./size(Mask, 1);

    NewMemory = [LastState,action,reward,CurrentState];
    Memory    = [Memory;NewMemory];

    if size(Memory,1) > 1000
        Memory = Memory(end - 500,:);
    end 
end
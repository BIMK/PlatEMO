function Memory = UpdateMemory(Memory,action,LastPopulation,LastMask,Population,Mask)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    LastPopulationHV = HV(LastPopulation,[3,3]);
    PopulationHV     = HV(Population,[3,3]);
    reward           = (PopulationHV - LastPopulationHV)/LastPopulationHV;
    LastSparse       = mean(LastMask(:));
    LastStd          = std(sum(LastMask,2));
    LastSparseDistuibution = histcounts(sum(LastMask,2)./size(LastMask,2),0:0.1:1)./size(LastMask,2);
    LastState        = [LastSparse,LastStd,LastSparseDistuibution];

    CurrentSparse = mean(Mask(:));
    CurrentStd    = std(sum(LastMask,2));
    CurrentSparseDistuibution = histcounts(sum(Mask,2)./size(Mask,2),0:0.1:1)./size(Mask,2);
    CurrentState  = [CurrentSparse,CurrentStd,CurrentSparseDistuibution];
    
    NewMemory = [LastState,action,reward,CurrentState];
    Memory    = [Memory;NewMemory];
end
function [Index,MAX] = CpuGroup(numberOfGroups,xPrime,numberOfVariables)       

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    varsPerGroup = floor(numberOfVariables/numberOfGroups);
    if varsPerGroup == 1
        Index = linspace(1,numberOfVariables,numberOfVariables);
        MAX   = numberOfVariables;
    else
        B      = ones(1,varsPerGroup*numberOfGroups);
        remain = ones(1,(numberOfVariables-varsPerGroup*numberOfGroups))*(numberOfGroups+1);
        R      = reshape(B,varsPerGroup,numberOfGroups);
        k      = linspace(1,numberOfGroups,numberOfGroups);
        index  = R.*repmat(k,varsPerGroup,1);
        index  = reshape(index,1,varsPerGroup*numberOfGroups);
        INDEX  = [index remain];
        [~,I]  = sort(xPrime);
        Index(I) = INDEX;
        if(mod(numberOfVariables,numberOfGroups)==0)
            MAX = numberOfGroups;
        else
            MAX = numberOfGroups + 1;
        end
    end
end
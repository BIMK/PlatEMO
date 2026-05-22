function [outIndexArray,numberOfGroupsArray] = CreateGroups(numberOfGroups, xPrime,numberOfVariables, method)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    outIndexArray       = [];
    numberOfGroupsArray = [];
    [noOfSolutions,D]   = size(xPrime);
    for sol = 1 : noOfSolutions
        switch method
            case 1 % linear grouping
                varsPerGroup = floor(numberOfVariables/numberOfGroups);
                outIndexList = [];
                for i = 1 : numberOfGroups-1
                   outIndexList = [outIndexList, ones(1,varsPerGroup).*i];
                end
                outIndexList = [outIndexList, ones(1,numberOfVariables-size(outIndexList,2)).*numberOfGroups];
            case 2 % orderByValueGrouping
                varsPerGroup = floor(numberOfVariables/numberOfGroups);
                vars         = xPrime(sol,:);
                [~,I]        = sort(vars);
                outIndexList = ones(1,numberOfVariables);
                for i = 1 : numberOfGroups-1
                   outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
                end
                outIndexList(I(((numberOfGroups-1)*varsPerGroup)+1:end)) = numberOfGroups;
            case 3 % random Grouping
                varsPerGroup = floor(numberOfVariables/numberOfGroups);
                outIndexList = [];
                for i = 1 : numberOfGroups-1
                   outIndexList = [outIndexList, ones(1,varsPerGroup).*i];
                end
                outIndexList = [outIndexList, ones(1,numberOfVariables-size(outIndexList,2)).*numberOfGroups];
                outIndexList = outIndexList(randperm(length(outIndexList)));
        end
        outIndexArray       = [outIndexArray;outIndexList];
        numberOfGroupsArray = [numberOfGroupsArray;numberOfGroups];
    end
end
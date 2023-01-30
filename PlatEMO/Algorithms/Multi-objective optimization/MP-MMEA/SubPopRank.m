function [FNmean,FNbest] = SubPopRank(Populations)
% Calculate the mean and best front number of each subpopulation

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Population = [Populations{:}];
    K          = length(Populations);
    Flag       = [];
    j          = 1;
    for i = 1 : K
        Flag(j:j+length(Populations{i})-1) = i;
        j = j + length(Populations{i});
    end
    FrontNoAll = NDSort(Population.objs,inf);
    FNmean     = zeros(1,K);
    FNbest     = zeros(1,K);
    for i = 1 : K
        FNmean(i) = mean(FrontNoAll(Flag==i));
        FNbest(i) = min(FrontNoAll(Flag==i));
    end
end
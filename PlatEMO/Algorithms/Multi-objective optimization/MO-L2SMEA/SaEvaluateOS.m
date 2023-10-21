function [objs, stds] = SaEvaluateOS(subproblemList,decs)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,nLinear] = size(decs);
    objs = zeros(N,nLinear);
    stds = zeros(N,nLinear);
    for k = 1 : nLinear
        [objs(:,k), stds(:,k)] = subproblemList{k}.fobj(decs(:,k));
    end
end
function [Population,LowDecs,MapCon,MapDiv,LB,UB] = ReMO_AdjustMap(Problem,HighDecs,ConVars,DivVars,ConDim,DivDim)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    MapCon = randn(ConDim,numel(ConVars));
    MapDiv = randn(DivDim,numel(DivVars));

    % Boundary
    LB = -1*ones(1,ConDim+DivDim);
    UB = 1*ones(1,ConDim+DivDim);

    LowDecs    = rand(Problem.N,ConDim+DivDim).*repmat(UB-LB,Problem.N,1)+repmat(LB,Problem.N,1);
    Population = ReMO_MapBack(Problem,LowDecs,MapCon,MapDiv,ConVars,DivVars);
end
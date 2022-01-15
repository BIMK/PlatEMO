function [df,dg] = FiniteDifference(X)
% Estimate the gradient of objective and constraints by finite difference

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    X1 = SOLUTION(repmat(X.dec,length(X.dec),1)+eye(length(X.dec))*1e-6);
    df = (X1.objs-repmat(X.obj,length(X1),1))/1e-6;
    dg = (X1.cons-repmat(X.con,length(X1),1))/1e-6;
end
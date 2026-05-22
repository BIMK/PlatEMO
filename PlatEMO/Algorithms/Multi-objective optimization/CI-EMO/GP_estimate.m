function [PopObj,PopObj_b,MSE,MSE_b] = GP_estimate(X,Model,numOBJ,beta)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

    [N,~]  = size(X);
    PopObj = zeros(N,numOBJ);
    MSE    = zeros(N,numOBJ);
    for i = 1: N
        for j = 1 : numOBJ
            [PopObj(i,j),~,MSE(i,j)] = predictor(X(i,:),Model{j});
        end
    end
    MSE = max(MSE,0);
    S_  = sqrt(MSE);

    MSE      = S_; % AUCB
    PopObj_b = PopObj;
    MSE_b    = MSE;
    PopObj   = PopObj+beta*MSE; % AUCB
end
function [W,W0] = setWeight(W,p)
% Set the weight vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    if 0 <= p && p <= 1
        %% Static approach
        W0 = [];
        M = size(W,2);
        % Distribution control of weight vector set
        TF = W < 1.0 / M;
        W(TF) = W(TF) * p * M;
        W(~TF) = 1.0 - (1.0 - W(~TF)) * (1.0 - p) * M / (M - 1);
    else
        %% Dynamic approach
        W0 = W;
    end
end 
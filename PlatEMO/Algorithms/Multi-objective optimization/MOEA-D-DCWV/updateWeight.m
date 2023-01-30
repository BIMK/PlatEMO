function W = updateWeight(objs,W)
% Update the weight vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    % Calculate the intermediate objective value p
    M       = size(objs,2);
    objs    = normalize(objs,'range');
    normP   = sqrt(sum(objs.^2,2));
    CosineP = sum(objs./M,2).*sqrt(M)./normP;
    [~,I]   = min(normP.*sqrt(1-CosineP.^2));
    p       = normP(I)*CosineP(I) / sqrt(M);
    % Distribution control of weight vector set
    TF = W < 1.0 / M;
    W(TF) = W(TF) * p * M;
    W(~TF) = 1.0 - (1.0 - W(~TF)) * (1.0 - p) * M / (M - 1);
end 
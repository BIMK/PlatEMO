function objhat = estimatePF(obj,validator,estimator)
% Estimate the Pareto front

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    %% X as input, Y as output
    obj = unique(normalize(obj,'range'),'rows');
    Y = vecnorm(obj,1,2);
    X = obj./Y;

    %% Generate the input W
    W = UniformPoint(2e4,width(obj),'ILD');
    W = validator(W,X);

    %% Estimate the objhat
    objhat = estimator(W,X,Y);
    objhat = objhat(NDSort(objhat,1)==1,:);
end 
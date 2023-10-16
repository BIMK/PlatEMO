function [theta,preApp,preCrowd] = UpdateTheta(preApp,preCrowd,newApp,newCrowd)
% Adaptive update function

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    appNew   = mean(newApp);
    crowdNew = mean(newCrowd);
    % calculate the rates of change of converegnce and diversity
    rateApp   = abs(appNew - preApp)/abs(preApp);
    rateCrowd = abs(crowdNew - preCrowd)/abs(preCrowd);
    % update theta
    Z          = rateApp - rateCrowd;
    theta      = 1 - exp(-exp(Z));
    preApp     = appNew;
    preCrowd   = crowdNew;
end
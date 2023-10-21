function refPoints = SelectRef(nLinear,xmean,sigma,B,Diag)
% This function is used to select several pair of reference points

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    D = length(xmean);
    % Select reference points
    refPoints = zeros(2*nLinear,D);
    for k = 1 : 2*nLinear
        refPoints(k,:) = (xmean + sigma * B * (Diag.*randn(D,1)))';
    end
    % Repair invalidate reference points, default BU and BD is 1,0
    refPoints(refPoints<0) = 1e-6;
    refPoints(refPoints>1) = 1-(1e-6);
end
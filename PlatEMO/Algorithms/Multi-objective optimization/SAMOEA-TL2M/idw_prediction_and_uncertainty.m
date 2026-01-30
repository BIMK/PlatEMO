function [prediction, s] = idw_prediction_and_uncertainty(x, known_points, values, p)
% IDW_PREDICTION_AND_UNCERTAINTY Calculate IDW prediction and uncertainty
%
%   [PREDICTION, S2] = IDW_PREDICTION_AND_UNCERTAINTY(X, KNOWN_POINTS,
%   VALUES, P) returns the IDW prediction and local variance estimate at
%   point X.
%
%   X is a vector representing the coordinates of the prediction point.
%   KNOWN_POINTS is an NxM matrix where each row represents the coordinates
%   of a known point.
%
%   VALUES is an Nx1 vector of values at the known points.
%   P is the power parameter for IDW, default is 2.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yuanchao Liu (email:liuyuanchao@ise.neu.edu.cn)
    
    if nargin < 4
        p = 2; % Default power parameter
    end
    % Number of known points
    n = size(x, 1);
    % Calculate the weights
    for i = 1 : n
        A          = repmat(x(i,:),size(known_points, 1),1) - known_points;
        distance   = sqrt(sum(A.^2, 2));
        weights    = 1 ./ (distance.^p);
        weights    = weights / sum(weights);
        prediction = sum(weights .* values);
        s(i)       = sum(weights .* (values - prediction).^2);
    end
end
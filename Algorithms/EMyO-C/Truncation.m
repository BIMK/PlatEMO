function [index] = Truncation(F, Z, remaining)
% Truncation of EMyO/C, based on clustering

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

% get size
[n, m] = size(F);

% calculate distances to reference point
Dist = sqrt(sum((F - repmat(Z, n, 1)).^2, 2));

% project to the hyperplane
F = (F - repmat(Z, n, 1))./repmat(sum(F - repmat(Z, n, 1), 2), 1, m);

% compute clusters
[C] = Clustering(F, remaining);

% find number of individuals in each cluster
numInd = zeros(numel(C), 1);
for i = 1 : numel(C)
    numInd(i) = numel(C{i});
end

index = zeros(remaining,1);
for i = 1 : remaining
    
    % get ind from cluster
    tempInd = C{i};
    
    % get individuals' distances to reference point 
    tempDist = Dist(tempInd);
    
    % find cluster's representative
    [~, idx] = min(tempDist);
    index(i) = tempInd(idx(1));
    
end

end
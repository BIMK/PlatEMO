function [C] = Clustering(data, k, varargin)
% Perform hierarchical agglomerative clustering to group points in
% matrix data into k clusters

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

% get size of data
n = size(data, 1);

% compute distances between points
distPoints = zeros(n);
for i = 1:n
    distPoints(i,:) = sqrt(sum(power(repmat(data(i,:),n,1)-data, 2), 2));
    distPoints(i,i) = inf;
end

% initialize distances between clusters
distClusters = distPoints;

% initially, each points belongs to a distinct cluster
C = cell(n, 1);
for i = 1:n;
    C{i} = i;
end

update = 0;
toRemove = zeros(n-k, 1);
for i = 1:n-k
    
    % update distances between clusters
    if update>0
        
        for jj = 1:n
            
            if isempty(C{jj}) || jj == update
                continue
            end
            
            % calculate distances between clusters
            distClusters(update, jj) = sum(sum(distPoints(C{update}, C{jj})))/(numel(C{update})*numel(C{jj}));
            
        end
        
    end
    
    % find two clusters with min distance
    [a, b] = find(distClusters==min(min(distClusters)));
    
    % merge these clusters
    C{a(1)} = horzcat(C{a(1)}, C{b(1)});
    
    % removed cluster
    toRemove(i) = b(1);
    C{toRemove(i)} = [];
    update = a(1);
    
    % remove from distance matrix
    distClusters(toRemove(i), :) = inf;
    distClusters(:, toRemove(i)) = inf;
    
end

% final clusters
C(toRemove) = [];
end
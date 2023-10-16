function varargout = cascade_cluster(P, Z, MtcStr, N, cat_flag)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% split the population
[O, ia, ~] = unique([P.objs], 'rows');
P = P(ia);
[index_frontier, NF_index] = find_frontiers(O, [P.cons]);
O = normalization(O);
%% attach frontiers to nearest reference vectors
if strcmp(MtcStr, 'PBI')
    error_type = 'distance';
elseif strcmp(MtcStr, 'PDM')
    error_type = 'sin';
end
[min_metric, allocation] = pair(O(index_frontier, :), Z, error_type);
%% identify the active reference vectors and create clusters
index_active_cluster = unique(allocation);
NC = numel(index_active_cluster);
cluster.center = []; cluster.queue_select = [];
clusters = repmat(cluster, NC, 1);
%% calculate metric and sort the frontier solutions
for i = 1: NC
    I = index_active_cluster(i);
    small_index = (allocation == I);
    index_cluster_frontier = index_frontier(small_index);
    if numel(index_cluster_frontier) == 1
        index_ascend = 1;
    else
        if strcmp(MtcStr, 'PBI')
            F_metric = 5 * min_metric(small_index)' + sum(O(index_cluster_frontier, :) .* repmat(Z(I, :), numel(index_cluster_frontier), 1), 2)' ./ repmat(norm(Z(I, :)), 1, numel(index_cluster_frontier));
        elseif strcmp(MtcStr, 'PDM')
            F_metric = 5 * (min_metric(small_index) .* sqrt(sum(O(index_cluster_frontier, :) .^ 2, 2)))' + mean(O(index_cluster_frontier, :), 2)';
        end
        if NC >= N
            [~, index_ascend] = min(F_metric);
        else
            [~, index_ascend] = sort(F_metric, 'ascend');
        end
    end
    clusters(i).queue_select = index_cluster_frontier(index_ascend);
    clusters(i).center = clusters(i).queue_select(1);
end
%% Attach the non-frontiers and sort
center_index = [clusters.center];
if nargout == 2
    CF = P(center_index);
end
if NC > N
    P = crowding_pick(P(center_index), N, 'precise');
elseif NC == N
    P = P(center_index);
else
    if cat_flag || numel(index_frontier) < N
        C = O(center_index, :);
        Distance2Centers = pdist2(O(NF_index, :), C, 'euclidean');
        [min_distance, allocation] = min(Distance2Centers, [], 2);
        active_clusters_with_NF_index = unique(allocation);
        for i = 1: numel(active_clusters_with_NF_index)
            I = active_clusters_with_NF_index(i);
            Cluster_NF_index = NF_index(allocation == I);
            NF_distance = min_distance(allocation == I);
            [~, index_ascend] = sort(NF_distance, 'ascend');
            clusters(I).queue_select = [clusters(I).queue_select, Cluster_NF_index(index_ascend)];
        end
    end
    %% Round-robin Picking
    j = 1;
    next_pop_index = zeros(1, N);
    clusters = clusters(randperm(NC));
    for pointer = 1: N
        while isempty(clusters(j).queue_select)
            j = mod(j, NC) + 1;
        end
        next_pop_index(pointer) = clusters(j).queue_select(1);
        clusters(j).queue_select(1) = [];
        j = mod(j, NC) + 1;
    end
    P = P(next_pop_index);
end
varargout{1} = P;
if nargout == 2
    varargout{2} = CF;
elseif nargout == 3
    varargout{2} = index_active_cluster;
    varargout{3} = setdiff(1: size(Z, 1), index_active_cluster);
end
end
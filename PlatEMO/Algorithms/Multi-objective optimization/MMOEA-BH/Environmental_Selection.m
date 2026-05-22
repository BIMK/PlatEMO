function [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection(Method, Problem, Archive, Population, LBA, PBA, n_PBA, varargin)
% Environmental_Selection - Unified entry point for environmental selection strategies.
% Method: 'APC', 'DBSCAN', or 'kmeans'
% varargin: Additional parameters needed by specific methods (e.g., maxCluster, N, epsilon, minpts)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    switch upper(Method)
        case 'APC'
            [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_APC(Problem, Archive, Population, LBA, PBA, n_PBA);
        case 'DBSCAN'
            % DBSCAN requires: N, maxCluster, epsilon(optional), minpts(optional)
            N          = varargin{1};
            maxCluster = varargin{2};
            if length(varargin) >= 4
                epsilon = varargin{3};
                minpts  = varargin{4};
                [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_DBSCAN(Problem, Archive, Population, LBA, PBA, n_PBA, N, maxCluster, epsilon, minpts);
            else
                [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_DBSCAN(Problem, Archive, Population, LBA, PBA, n_PBA, N, maxCluster);
            end
        case 'KMEANS'
            % KMEANS requires: maxCluster
            maxCluster = varargin{1};
            [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_kmeans(Problem, Archive, Population, LBA, PBA, n_PBA, maxCluster);
        otherwise
            error('Unknown Environmental Selection Method: %s', Method);
    end
end

function [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_APC(Problem, Archive, Population, LBA, PBA, n_PBA)
    N         = Problem.N;
    temp_p    = [Population, Archive, LBA{:}];
    temp_decs = temp_p.decs;
    [~,idx]   = unique(temp_decs(:,1));
    Archive   = temp_p(idx);
    
    temp_Archive      = Archive;
    [idx, maxCluster] = MMOEA_Utils.APC(temp_Archive.decs);
    temp_Archive_nd   = [];
    temp_idx_d_nd     = [];
    
    for i = 1 : maxCluster
        indv_t1         = temp_Archive(idx==i);
        num_ND          = NDSort(indv_t1.objs,1)==1;
        temp_Archive_nd = [temp_Archive_nd,indv_t1(num_ND)];
        temp_idx_d_nd   = [temp_idx_d_nd;i*ones(sum(num_ND),1)];
    end
    Archive = temp_Archive_nd;
    idx     = temp_idx_d_nd;
    
    while numel(Archive)>N
        num_cluster = tabulate(idx);
        [~,idx_max] = max(num_cluster(:,2));
        density_y   = MMOEA_Utils.calc_pccs(Archive.objs);
        sel_cluster = find(idx == idx_max);
        density_yt  = density_y(sel_cluster);
        [~,idx1]    = sort(density_yt,'descend');
        try
            idx2 = idx1(1:Problem.M+1);
        catch
            idx2 = idx1;
        end
        density_x  = MMOEA_Utils.calc_pccs(Archive.decs);
        density_xt = density_x(sel_cluster(idx2));
        [~,idx3]   = max(density_xt);
        Archive(sel_cluster(idx2(idx3))) = [];
        idx(sel_cluster(idx2(idx3)))     = [];
    end
    
    a1    = tabulate(idx);
    idx_t = sort(find(a1(:,2)==0),'descend')';
    for i = idx_t
        idx(idx>i,:) = idx(idx>i,:)-1;
    end
    maxCluster = max(idx);
    
    for i = 1 : N
        [~, temp_PBA,temp_PBA_SCD] = MMOEA_Utils.nd_pccs_sort([Population(i), PBA{i}]);
        PBA{i}                     = temp_PBA(1:min(numel(temp_PBA),n_PBA));
        PBA_SCD{i}                 = temp_PBA_SCD(1:min(numel(temp_PBA),n_PBA));
    end
    for i = 1 : maxCluster
        [~, LBA{i}, LBA_SCD{i}] = MMOEA_Utils.nd_pccs_sort(Archive(i==idx),ceil(N/maxCluster));
    end
    LBA(maxCluster+1:end)     = [];
    LBA_SCD(maxCluster+1:end) = [];
end

function [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_DBSCAN(Problem, Archive, Population, LBA, PBA, n_PBA, N, maxCluster, epsilon, minpts)
    if ~exist('epsilon','var'), [epsilon, minpts] = deal(0.15, 5); end
    
    temp_p    = [Population, Archive, LBA{:}];
    temp_decs = temp_p.decs;
    [~,idx]   = unique(temp_decs(:,1));
    Archive   = temp_p(idx);
    
    temp_Archive = Archive;
    idx_1        = dbscan(normalize(temp_Archive.decs, 'range'), epsilon, minpts);
    if ismember(-1, idx_1), maxCluster_1 = numel(unique(idx_1))-1; else, maxCluster_1 = numel(unique(idx_1)); end
    
    idx_2        = dbscan(normalize(temp_Archive.objs, 'range'), epsilon, minpts);
    if ismember(-1, idx_2), maxCluster_2 = numel(unique(idx_2))-1; else, maxCluster_2 = numel(unique(idx_2)); end
    
    temp_num = tabulate(idx_1);
    if ismember(-1, idx_1), temp_num = temp_num(2:end,2); else, temp_num = temp_num(:,2); end
    if max(temp_num)/min(temp_num) > 5, flag_dec = 0; else, flag_dec = 1; end
    
    if maxCluster_1==1 || maxCluster_2==1
        if maxCluster_1==1, flag_dec_one = 1; else, flag_dec_one = 0; end
        if flag_dec_one, [idx, maxCluster_d] = deal(idx_2, maxCluster_2); else, [idx, maxCluster_d] = deal(idx_1, maxCluster_1); end
    else
        if maxCluster_1 >= maxCluster_2 && flag_dec, [idx, maxCluster_d] = deal(idx_1, maxCluster_1); else, [idx, maxCluster_d] = deal(idx_2, maxCluster_2); end
    end
    
    temp_Archive_nd = [];
    for i = -1 : maxCluster_d
        indv_t1         = temp_Archive(idx==i);
        temp_Archive_nd = [temp_Archive_nd,indv_t1(NDSort(indv_t1.objs,1)==1)];
    end
    Archive = temp_Archive_nd;
    
    idx_1 = dbscan(normalize(Archive.decs, 'range'), epsilon, minpts);
    idx_2 = dbscan(normalize(Archive.objs, 'range'), epsilon, minpts);
    if ismember(-1, idx_1), maxCluster_1 = numel(unique(idx_1))-1; else, maxCluster_1 = numel(unique(idx_1)); end
    if ismember(-1, idx_2), maxCluster_2 = numel(unique(idx_2))-1; else, maxCluster_2 = numel(unique(idx_2)); end
    
    if maxCluster_1 == maxCluster_2, idx_del = any([idx_1==-1,idx_2==-1],2); else, idx_del = all([idx_1==-1,idx_2==-1],2); end
    if sum(idx_del==1) == numel(idx_del), idx_del = []; end
    Archive(idx_del) = [];
    
    idx_1 = dbscan(normalize(Archive.decs, 'range'), epsilon, minpts);
    idx_2 = dbscan(normalize(Archive.objs, 'range'), epsilon, minpts);
    if maxCluster_1 >= maxCluster_2, idx = idx_1; else, idx = idx_2; end
    
    while numel(Archive) > N
        num_cluster = tabulate(idx);
        [~,idx_max] = max(num_cluster(:,2));
        idx_max     = num_cluster(idx_max,1);
        density_y   = MMOEA_Utils.calc_pccs(Archive.objs);
        sel_cluster = find(idx == idx_max);
        density_yt  = density_y(sel_cluster);
        [~,idx1]    = sort(density_yt,'descend');
        try
            idx2 = idx1(1:Problem.M+1);
        catch
            idx2 = idx1;
        end
        density_x  = MMOEA_Utils.calc_pccs(Archive.decs);
        density_xt = density_x(sel_cluster(idx2));
        [~,idx3]   = max(density_xt);
        Archive(sel_cluster(idx2(idx3))) = [];
        idx(sel_cluster(idx2(idx3)))     = [];
    end
    
    maxCluster = numel(unique(idx));
    if ismember(-1, idx)
        idx(idx==-1) = maxCluster;
    end
    
    for i = 1 : numel(Population)
        [~, temp_PBA,temp_PBA_SCD] = MMOEA_Utils.nd_pccs_sort([Population(i), PBA{i}]);
        PBA{i}                     = temp_PBA(1:min(numel(temp_PBA),n_PBA));
        PBA_SCD{i}                 = temp_PBA_SCD(1:min(numel(temp_PBA),n_PBA));
    end
    for i = 1 : maxCluster
        [~, LBA{i}, LBA_SCD{i}] = MMOEA_Utils.nd_pccs_sort(Archive(i==idx),ceil(N/maxCluster));
    end
    LBA(maxCluster+1:end)     = [];
    LBA_SCD(maxCluster+1:end) = [];
end

function [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_kmeans(Problem, Archive, Population, LBA, PBA, n_PBA, maxCluster)
    N         = Problem.N;
    temp_p    = [Population, Archive, LBA{:}];
    temp_decs = temp_p.decs;
    [~,idx_1] = unique(temp_decs(:,1));
    Archive   = temp_p(idx_1);
    
    temp_Archive    = Archive;
    idx             = kmeans(Archive.decs, maxCluster);    
    temp_Archive_nd = [];
    temp_idx_d_nd   = [];
    
    for i = 1 : maxCluster
        indv_t1         = temp_Archive(idx==i);
        num_ND          = NDSort(indv_t1.objs,1)==1;
        temp_Archive_nd = [temp_Archive_nd,indv_t1(num_ND)];
        temp_idx_d_nd   = [temp_idx_d_nd;i*ones(sum(num_ND),1)];
    end
    Archive = temp_Archive_nd;
    idx     = temp_idx_d_nd;
    
    while numel(Archive) > N
        num_cluster = tabulate(idx);
        [~,idx_max] = max(num_cluster(:,2));
        density_y   = MMOEA_Utils.calc_pccs(Archive.objs);
        sel_cluster = find(idx == idx_max);
        density_yt  = density_y(sel_cluster);
        [~,idx1]    = sort(density_yt,'descend');
        try
            idx2 = idx1(1:Problem.M+1);
        catch
            idx2 = idx1;
        end
        density_x  = MMOEA_Utils.calc_pccs(Archive.decs);
        density_xt = density_x(sel_cluster(idx2));
        [~,idx3]   = max(density_xt);
        Archive(sel_cluster(idx2(idx3))) = [];
        idx(sel_cluster(idx2(idx3)))     = [];
    end
    
    a1    = tabulate(idx);
    idx_t = sort(find(a1(:,2)==0),'descend')';
    for i = idx_t
        idx(idx>i,:) = idx(idx>i,:)-1;
    end
    maxCluster = max(idx);
    
    for i = 1 : N
        [~, temp_PBA,temp_PBA_SCD] = MMOEA_Utils.nd_pccs_sort([Population(i), PBA{i}]);
        PBA{i}                     = temp_PBA(1:min(numel(temp_PBA),n_PBA));
        PBA_SCD{i}                 = temp_PBA_SCD(1:min(numel(temp_PBA),n_PBA));
    end
    for i = 1 : maxCluster
        [~, LBA{i}, LBA_SCD{i}] = MMOEA_Utils.nd_pccs_sort(Archive(i==idx),ceil(N/maxCluster));
    end
    LBA(maxCluster+1:end)     = [];
    LBA_SCD(maxCluster+1:end) = [];
end
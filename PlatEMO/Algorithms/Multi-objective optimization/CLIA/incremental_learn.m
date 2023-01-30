function [Z, SVM] = incremental_learn(Z, active_clusters_index, index_inactive_clusters, frontiers, SVM, Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

global MAX_REFERENCE_NUM_FLAG MAX_ARCHIVE_SIZE disable_flag normalization_str delta learning_initial_flag;
if learning_initial_flag == true
    learning_initial_flag = false;
    MAX_ARCHIVE_SIZE = ceil(0.33 * Problem.M * Problem.N);
end
if disable_flag && numel(frontiers) > 0.9 * MAX_ARCHIVE_SIZE || numel(active_clusters_index) > 0.95 * Problem.N && size(Z, 1) == Problem.N
    return;
end
status = check_status(size(Z, 1), active_clusters_index, index_inactive_clusters, Problem);
if strcmp(status, 'unstable') || numel(frontiers) < 0.9 * MAX_ARCHIVE_SIZE
    return;
end
if numel(active_clusters_index) > Problem.N && numel(active_clusters_index) == size(Z, 1)
    Z = generate_reference([], 0, Problem, -1, 'normal');
elseif numel(active_clusters_index) < 0.95 * Problem.N
    if MAX_REFERENCE_NUM_FLAG
        return;
    end
    [~, allocation] = pair([frontiers.objs], Z, 'sin');
    Z_active_old = Z(unique(allocation), :);
    [Z, ~] = generate_reference([frontiers.objs], numel(active_clusters_index), Problem, 1, 'normal');
    if strcmp(normalization_str, 'normalize')
        if size(Z, 1) > 4 * Problem.N
            normalization_str = '';
            [~, allocation] = pair([frontiers.objs], Z, 'sin');
            index_inactive_clusters = setdiff(1: size(Z, 1), unique(allocation));
            P = truncate([frontiers.objs]); N = truncate(Z(index_inactive_clusters, :));
            SVM = learn(P, N, SVM, 'svmtrain2');
            Z_next = reduce_useless(Z, delta, SVM, Problem);
        else
            Z_next = Z;
        end
    elseif size(Z, 1) > 2 * Problem.N
        picked_frontiers = crowding_pick(frontiers, 2 *  Problem.N, 'precise');
        [~, allocation] = pair([frontiers.objs], Z, 'sin');
        index_inactive_clusters = setdiff(1: size(Z, 1), unique(allocation));
        P = truncate([picked_frontiers.objs]); N = truncate(Z(index_inactive_clusters, :));
        if size(P, 1) > 0 && size(N, 1) > 0
            [~, y_active] = predict(P, SVM); p_index = find(y_active == 1);
            [~, y_inactive] = predict(N, SVM); n_index = find(y_inactive == -1);
            P = P(setdiff(1: size(P, 1), p_index), :); N = N(setdiff(1: size(N, 1), n_index), :);
            SVM = learn(P, N, SVM, 'svmtrain2');
            Z_next = reduce_useless(Z, delta, SVM, Problem);
        else
            Z_next = Z;
        end
    else
        Z_next = Z;
    end
    [~, allocation] = pair([frontiers.objs], Z, 'sin');
    Z = unique([Z_active_old; Z_next; Z(unique(allocation), :)], 'rows');
end
end

function Z = reduce_useless(Z, threshold, SVM, Problem)
original_size = size(Z, 1);
if original_size <= Problem.N
    index_remain = 1: original_size;
else
    [score, ~] = predict(Z, SVM);
    if length(score) == 1
        index_remain = 1: original_size;
    elseif isempty(threshold)
        index_remain = [];
    elseif threshold > 0 && threshold <= 1
        threshold = min(threshold, 0.5);
        index_remain = find(score >= threshold);
    elseif threshold > 1
        S = sort(score, 'descend');
        threshold = S(min(max(Problem.N, threshold), numel(S)));
        index_remain = find(score >= threshold);
    end
    if numel(index_remain) < Problem.N
        S = sort(score, 'descend');
        threshold = S(min(Problem.N, numel(S)));
        index_remain = find(score >= threshold);
    end
end
Z = Z(index_remain, :);
fprintf('reference vectors reduced from %d to %d (-%d)\n', original_size, numel(index_remain), original_size - numel(index_remain));
end

function [status, index_i_active, index_i_inactive] = check_status(total_number, active_index, inactive_index, Problem)
global normalization_str current_density status_initial_flag stable_threshold;
persistent density history consecutive_stable_counter threshold;
if status_initial_flag == true
    status_initial_flag = false;
    consecutive_stable_counter = 0;
    density = current_density;
    history = zeros(1, total_number);
    if Problem.M == 5
        pointer = 1;
    elseif Problem.M == 10
        pointer = 2;
    elseif Problem.M == 15
        pointer = 3;
    else
        pointer = 0;
    end
    if pointer > 0 && stable_threshold(pointer) > 0
        threshold = stable_threshold(pointer);
    else
        threshold = min(20, max(5, ceil(Problem.maxFE / 2e4)));
    end
end
if strcmp(normalization_str, 'normalize')
    fluctuation = 3;
else
    fluctuation = 1e-2;
end
if density ~= current_density
    density = current_density;
    status_initial_flag = false;
    consecutive_stable_counter = 0;
    history = zeros(1, total_number);
end
current = zeros(size(history));
current(active_index) = 1;
if norm(history - current) <= fluctuation
    consecutive_stable_counter = consecutive_stable_counter + 1;
else
    history = current;
    consecutive_stable_counter = 0;
end
if consecutive_stable_counter >= threshold
    consecutive_stable_counter = 0;
    threshold = max(5, threshold - 1);
    status = 'stable';
    index_i_active = active_index; index_i_inactive = inactive_index;
else
    status = 'unstable';
    index_i_active = []; index_i_inactive = [];
end
end

function SVM = learn(P, N, SVM, training_mode)
length_P = size(P, 1); length_N = size(N, 1);
if length_P == 0 || length_N == 0
    return;
end
X = project2simplex([P; N]); y = [ones(length_P, 1); -1 * ones(length_N, 1)];
if ~isempty(SVM.y_mer) % NOT NEW
    [SVM.a, SVM.b, SVM.g, SVM.ind, SVM.uind, SVM.X_mer, SVM.y_mer, SVM.Rs, SVM.Q] = feval(training_mode, X', y, SVM.C);
else
    [SVM.a, SVM.b, SVM.g, SVM.ind, SVM.uind, SVM.X_mer, SVM.y_mer, SVM.Rs, SVM.Q] = feval(training_mode, X', y, SVM.C, SVM.type, SVM.scale);
end
end

function projected_points = project2simplex(O)
global hp_initial_flag;
persistent V Q;
if isempty(V) || hp_initial_flag == true
    M = size(O, 2);
    Q = orth(eye(M, M) - repmat(ones(M, 1) ./ M, 1, M));
    hp_initial_flag = false;
end
if size(O, 1)
    projected_points = O * Q;
else
    projected_points = [];
end
end

function [score, class] = predict(X, SVM)
class = nan(1, size(X, 1));
score = class';
if ~isempty(SVM.y_mer)
    block_size = 10000;
    X = project2simplex(X);
    blocks = floor(size(X, 1) / block_size);
    score = nan(size(size(X, 1), 1));
    B = mat2cell(X, [block_size * ones(1, blocks), size(X, 1) - block_size * blocks], size(X, 2));
    for i = 1: blocks
        score((i - 1) * block_size + 1: i * block_size) = 0.5 + 0.5 * elliotsig(svmeval((B{i})', SVM.a, SVM.b, SVM.ind, SVM.X_mer, SVM.y_mer, SVM.type, SVM.scale));
    end
    score(block_size * blocks + 1: size(X, 1)) = 0.5 + 0.5 * elliotsig(svmeval((B{blocks + 1})', SVM.a, SVM.b, SVM.ind, SVM.X_mer, SVM.y_mer, SVM.type, SVM.scale));
    class(score >= 0.5) = 1; class(score < 0.5) = -1;
end
end

function Z = truncate(Z)
Z = Z ./ repmat(sum(Z, 2), 1, size(Z, 2));
end
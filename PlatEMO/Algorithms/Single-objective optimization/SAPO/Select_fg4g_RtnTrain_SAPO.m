function [pop, fit, con, train_x, train_y, train_c] = Select_fg4g_RtnTrain_SAPO(db_x, db_y, db_c, ind_infea_better, ind_worse, ncon, N, n, change_data_n, ch_d_n_thres, data_n, data_n_2, data_times)
%% Select solutions with good f and gm to improve gn

    i = 1;
    ind_pop_cell = cell(ncon, 1);
    ind_f = [];
    % f <= bf and gn > 0 and gm ascending order
    ind  = db_c(ind_infea_better, n) > 0;
    ind2 = ind_infea_better(ind);
    [~, ind3] = sort(db_c(ind2, n));
    ind_f = [ind_f; ind2(ind3)];
    
    % f > bf and gn > 0 and f ascending order
    ind  = db_c(ind_worse, n) > 0;
    ind2 = ind_worse(ind);
    [~, ind3] = sort(db_y(ind2, 1));
    ind_f = [ind_f; ind2(ind3)];

    ind_pop_cell{i} = ind_f.';

    for m = 1:ncon
        if m == n
            continue
        end
        i = i + 1;
        ind_g = [];
        % f <= bf and gm <= 0 and gn > 0 and gn ascending order
        ind  = db_c(ind_infea_better, m) <= 0;
        ind2 = ind_infea_better(ind);
        ind3 = db_c(ind2, n) > 0;
        ind4 = ind2(ind3);
        [~, ind5] = sort(db_c(ind4, n));
        ind_g = [ind_g; ind4(ind5)];
        
        % f > bf and gm <= 0 and gn > 0 and f ascending order
        ind  = db_c(ind_worse, m) <= 0;
        ind2 = ind_worse(ind);
        ind3 = db_c(ind2, n) > 0;
        ind4 = ind2(ind3);
        [~, ind5] = sort(db_y(ind4, 1));
        ind_g = [ind_g; ind4(ind5)];
        
        % gm > 0 and gn > 0 and gm ascending order
        ind2 = find(db_c(:, m) > 0);
        ind3 = db_c(ind2, n) > 0;
        ind4 = ind2(ind3);
        [~, ind5] = sort(db_c(ind4, m));
        ind_g = [ind_g; ind4(ind5)];
    
        ind_pop_cell{i} = ind_g.';
    end

    max_length = max(cellfun(@numel, ind_pop_cell));
    ind_pop_cell = cellfun(@(x) [x, zeros(1, max_length - numel(x))], ind_pop_cell, 'UniformOutput', false);
    ind_pop = vertcat(ind_pop_cell{:});
    ind_pop = ind_pop(:);
    ind_pop(ind_pop == 0) = [];
    ind_pop = unique(ind_pop, 'stable');

    ind_train = ind_pop;

    if size(ind_pop, 1) >= N
        ind_pop = ind_pop(1:N, :);
    end
    
    pop = db_x(ind_pop, :); fit = db_y(ind_pop, :); con = db_c(ind_pop, :); % v = max(0, con); cv = sum(v, 2);

    switch change_data_n
        case 0
            train_size = data_n;
        case 1
            train_size = data_n;
            if size(db_x, 2) >= ch_d_n_thres
                train_size = data_n_2;
            end
        case 2
            train_size = data_times * size(db_x, 2);
    end
    if size(ind_train, 1) >= train_size
        ind_train = ind_train(1:train_size, :);
    end
    train_x = db_x(ind_train, :); train_y = db_y(ind_train, :); train_c = db_c(ind_train, :); % train_v = max(0, train_c); train_cv = sum(train_v, 2);
end
function [pop, fit, cv, train_x, train_y, train_c] = Select_g4f_half_RtnTrain_SAPO(db_x, db_y, db_c, ind_infea_better, ind_worse, ncon, N, db_feasible, len, change_data_n, ch_d_n_thres, data_n, data_n_2, data_times)
%% Select solutions with good g to improve f

    ind_pop_cell = cell(ncon, 1);
    ind_pop_cell_train = cell(ncon, 1);
    for m = 1:ncon
        ind_g = [];
        % cv = 0 and f ascending order
        ind = find(db_feasible((1:len).'));
        [~, ind_tmp] = sort(db_y(db_feasible((1:len).'), 1));
        ind_g = [ind_g; ind(ind_tmp)];
        ind_g_train = ind_g;
        if size(ind_g, 1) >= floor(N/2)
            ind_g = ind_g(1:floor(N/2), :);
        end
        
        % f <= bf and gm <= 0 and gm ascending order
        ind  = db_c(ind_infea_better, m) <= 0;
        ind2 = ind_infea_better(ind);
        [~, ind3] = sort(db_c(ind2, m));
        ind_g = [ind_g; ind2(ind3)];
        ind_g_train = [ind_g_train; ind2(ind3)];
        
        % f > bf and gm <= 0 and f ascending order
        ind  = db_c(ind_worse, m) <= 0;
        ind2 = ind_worse(ind);
        [~, ind3] = sort(db_y(ind2, 1));
        ind_g = [ind_g; ind2(ind3)];
        ind_g = unique(ind_g, 'stable');
        ind_g_train = [ind_g_train; ind2(ind3)];
        ind_g_train = unique(ind_g_train, 'stable');
    
        % f <= bf and gm > 0 and gm ascending order
        ind  = db_c(ind_infea_better, m) > 0;
        ind2 = ind_infea_better(ind);
        [~, ind3] = sort(db_c(ind2, m));
        ind_g = [ind_g; ind2(ind3)];
        ind_g = unique(ind_g, 'stable');
        ind_g_train = [ind_g_train; ind2(ind3)];
        ind_g_train = unique(ind_g_train, 'stable');
        
        % f > bf and gm > 0 and f ascending order
        ind  = db_c(ind_worse, m) > 0;
        ind2 = ind_worse(ind);
        [~, ind3] = sort(db_y(ind2, 1));
        ind_g = [ind_g; ind2(ind3)];
        ind_g = unique(ind_g, 'stable');
        ind_g_train = [ind_g_train; ind2(ind3)];
        ind_g_train = unique(ind_g_train, 'stable');
        
        ind_pop_cell{m} = ind_g.';
        ind_pop_cell_train{m} = ind_g_train.';
    end

    max_length = max(cellfun(@numel, ind_pop_cell));
    ind_pop_cell = cellfun(@(x) [x, zeros(1, max_length - numel(x))], ind_pop_cell, 'UniformOutput', false);
    ind_pop = vertcat(ind_pop_cell{:});
    ind_pop = ind_pop(:);
    ind_pop(ind_pop == 0) = [];
    ind_pop = unique(ind_pop, 'stable');

    if size(ind_pop, 1) >= N
        ind_pop = ind_pop(1:N, :);
    end

    pop = db_x(ind_pop, :); fit = db_y(ind_pop, :); con = db_c(ind_pop, :); v = max(0, con); cv = sum(v, 2);

    max_length = max(cellfun(@numel, ind_pop_cell_train));
    ind_pop_cell_train = cellfun(@(x) [x, zeros(1, max_length - numel(x))], ind_pop_cell_train, 'UniformOutput', false);
    ind_train = vertcat(ind_pop_cell_train{:});
    ind_train = ind_train(:);
    ind_train(ind_train == 0) = [];
    ind_train = unique(ind_train, 'stable');
    
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
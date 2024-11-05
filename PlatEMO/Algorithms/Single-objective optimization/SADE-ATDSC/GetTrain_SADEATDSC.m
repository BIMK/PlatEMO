function [train_x, train_y] = GetTrain_SADEATDSC(id, db_x, db_y, pop, N, data_n)
    switch id
        case 1  % All Data
            train_x = db_x; train_y = db_y;
        case 2  % Top
            [B, I]  = sort(db_y); 
            train_x = db_x(I(1:N), :); train_y = B(1:N);
        case 3  % Recent
            train_x = db_x(end-N+1:end, :); train_y = db_y(end-N+1:end); 
        case 4  % Neighbor
            dists   = pdist2(pop, db_x); [~, I] = sort(dists, 2); I = unique(I(:, 1:min(data_n, size(db_x, 1))));
            train_x = db_x(I, :); train_y = db_y(I);
        case 5  % Near Best
            [~, I]  = sort(db_y); b1 = db_x(I(1), :); dist = pdist2(b1, db_x); [~, I] = sort(dist);
            train_x = db_x(I(1:N), :); train_y = db_y(I(1:N));
    end
end
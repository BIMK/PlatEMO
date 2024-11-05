function [bx, by, bc, bcv, bf] = Findbest_SAPO(db_x, db_y, db_c, db_cv)
%% Get the best feasible solution

    ind = find(db_cv == 0);
    if ~isempty(ind)    % feasible
        [bf, bid] = min(db_y(ind, 1));
        bx  = db_x(ind(bid), :);
        by  = db_y(ind(bid), :);
        bc  = db_c(ind(bid), :);
        bcv = 0;
    else                % infeasible
        [~, bid] = min(db_cv);
        bx  = db_x(bid,  :);
        by  = db_y(bid,  :);
        bc  = db_c(bid,  :);
        bcv = db_cv(bid, :);
        bf  = NaN;
    end
end
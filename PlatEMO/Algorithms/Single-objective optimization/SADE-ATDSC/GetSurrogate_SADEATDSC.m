function surr = GetSurrogate_SADEATDSC(train_x, train_y, RBFN_lib, spread, kernel)
    switch RBFN_lib
        case 1  % newrbe
            if isnan(spread)
                D    = size(train_x, 2);
                pair = pdist2(train_x, train_x);
                spr  = max(max(pair,[],2)) * (D * size(train_x, 1)) .^ (-1 / D);
            else
                spr = spread;
            end
            surr = newrbe(train_x.', train_y.', spr);
        case 2  % RBF + RBF_eval
            surr = struct; surr.train_x = train_x;
            [surr.lambda, surr.gamma, surr.rho] = RBF(train_x, train_y, kernel);
        case 3  % rbfcreate + rbfinterp
            surr = rbfcreate(train_x.', train_y.', 'RBFFunction', kernel);
    end
end
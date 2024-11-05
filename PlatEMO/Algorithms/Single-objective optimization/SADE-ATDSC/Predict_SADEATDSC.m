function pred_Fs = Predict_SADEATDSC(x, surr, RBFN_lib, kernel)
    switch RBFN_lib
        case 1  % newrbe
            pred_Fs = sim(surr, x.').';
        case 2  % RBF + RBF_eval
            pred_Fs = RBF_eval(x, surr.train_x, surr.lambda, surr.gamma, surr.rho, kernel);
        case 3  % rbfcreate + rbfinterp
            pred_Fs = rbfinterp(x.', surr).';
    end
end
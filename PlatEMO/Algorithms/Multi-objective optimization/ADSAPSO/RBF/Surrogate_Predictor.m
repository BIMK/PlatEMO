function Surrogate_Obj=Surrogate_Predictor(X, Surrogate_model, M)

% This function is written by Jianqing Lin

    %% Prediction of RBF
    N             = size(X,1);
    Surrogate_Obj = zeros(N,M);


    for j=1:M
        Surrogate_Obj(:,j) = RBFInterp(X, Surrogate_model{j});
    end
end

function y = RBFInterp(x, para)
    ax = para.nodes;
    nx = size(x, 1);
    np = size(ax, 1);    % np: the size of data set

    xmin = para.xmin;
    xmax = para.xmax;
    ymin = para.ymin;
    ymax = para.ymax;

    % normalization
    x = 2./(repmat(xmax - xmin, nx, 1)) .* (x - repmat(xmin, nx, 1)) - 1;

    r = dist(x, ax');
    switch para.kernel
        case 'gaussian'
    %         Phi = 1 / sqrt(2 * pi) * exp(- 1/2 * r.^2);
            Phi = radbas(sqrt(-log(.5))*r);
        case 'cubic'
            Phi = r.^3;
    end

    y = Phi * para.alpha + [ones(nx, 1), x] * para.beta;
    % renormalization
    y = repmat(ymax - ymin, nx, 1)./2 .* (y + 1) + repmat(ymin, nx, 1);
end
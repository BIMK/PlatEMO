function y = RBFInterp_SAPO(x, para)
    ax = para.nodes;
    nx = size(x, 1);
    
    xmin = para.xmin;
    xmax = para.xmax;
    ymin = para.ymin;
    ymax = para.ymax;
    % normalization
    for i = 1 : length(xmin)
       if xmin(i) ~= xmax(i)
           x(:, i) = 2 ./ (repmat(xmax(i) - xmin(i), nx, 1)) .* (x(:, i) - repmat(xmin(i), nx, 1)) - 1;
       end
    end
    
    r = dist(x, ax');
    switch para.kernel
        case 'gaussian'
            Phi = radbas(sqrt(-log(.5))*r);
        case 'cubic'
            Phi = r.^3;
        case 'multiquadric'
            Phi=sqrt(r.^2+0.8^2);
    end
    
    y = Phi * para.alpha + [ones(nx, 1), x] * para.beta;
    % renormalization
    for i = 1 : length(ymin)
       if ymin(i) ~= ymax(i)
           y(:, i) = repmat(ymax(i) - ymin(i), nx, 1)./2 .* (y(:, i) + 1) + repmat(ymin(i), nx, 1);
       end
    end
end
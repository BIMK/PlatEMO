function f = rbfinterp(x, options)

    phi      = options.('rbfphi');
    rbfconst = options.('RBFConstant');
    nodes    = options.('x');
    rbfcoeff = (options.('rbfcoeff'))';
    
    [dim              n] = size(nodes);
    [dimPoints  nPoints] = size(x);
    
    if (dim~=dimPoints)
        error(sprintf('x should have the same number of rows as an array used to create RBF interpolation'));
    end
    
    f = zeros(1, nPoints);
    
    for i = 1 : nPoints
        s = 0;
        r =  (x(:,i)*ones(1,n)) - nodes;
        r = sqrt(sum(r.*r, 1));
        s = rbfcoeff(n+1) + sum(rbfcoeff(1:n).*feval(phi, r, rbfconst));
        for k = 1 : dim
            s = s+rbfcoeff(k+n+1)*x(k,i);     % linear part
        end
        f(i) = s;
    end
end
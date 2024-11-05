function [Y,Error] = Sammon(X,d)
% Sammon mapping on X
% Usage: [Y,Error] = sammon(X,D)

% References:
% Sammon, John W. Jr., "A Nolinear Mapping for Data Structure Analysis",
% IEEE Transactions on Computers, vol. C-18, no. 5, pp 401-409, May 1969.

    %% Default setting
    maxHalves = 20;
    maxIter   = 500;
    theta     = 1e-9;
    
    %% Obtain distance matrix
    Dist = pdist2(X,X);
    
    %% Remaining inialisation
    N     = size(X,1);
    scale = 0.5/sum(Dist(:));
    Dist  = Dist + eye(N);
    Dinv  = 1./Dist;
    
    %% Initialized by PCA, map to d dimension
    [UU,DD] = svd(X);
    Y       = UU(:,1:d)*DD(1:d,1:d);
    %% Randomly initialize Y
    % Y      = randn(N,d);
    
    % Obtain error
    oneMatrix = ones(N,d);
    dist      = pdist2(Y,Y) + eye(N);
    dinv      = 1./dist;
    delta     = Dist - dist;
    Error     = sum(sum((delta.^2).*Dinv));
    
    %% Optimize process
    for i = 1 : maxIter
        % Compute gradient, Hessian and search direction
        delta    = dinv - Dinv;
        deltaOne = delta*oneMatrix;
        g        = delta*Y - Y.*deltaOne;
        dinv3    = dinv.^3;
        Y2       = Y.^2;
        H        = dinv3*Y2 - deltaOne - 2*Y.*(dinv3*Y) + Y2.*(dinv3*oneMatrix);
        s        = -g(:)./abs(H(:));
        yOld     = Y;
        
        % Use step-halving procedure to ensure progress id made
        for j = 1 : maxHalves
            Y(:)   = yOld(:) + s;
            dist   = pdist2(Y,Y) + eye(N);
            dinv   = 1./dist;
            delta  = Dist - dist;
            newError = sum(sum((delta.^2).*Dinv));
            if newError < Error
                break;
            else
                s = 0.5*s;
            end
        end
        
        % Evaluate termination criterion
        if abs((Error-newError)/Error) < theta
            break;
        end
        
        % Update error
        Error = newError;
    end
    Error = Error * scale;
end
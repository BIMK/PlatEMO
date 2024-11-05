function Yest=RBF_eval(X,S,lambda,gamma,rho,flag)   % For Gaussian basis function
%--------------------------------------------------------------------------
%Copyright (c) 2012 by Juliane Mueller
%
% This file is part of the surrogate model module toolbox.
%
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%Tampere University of Technology, Finland
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%input:
%X are points where function values should be calculated
%S are points where the function values are known
%lambda parameter vector
%gamma contains the parameters of the optional polynomial tail
%flag is a string indicating which RBF to be used
%output: 
%the estimated function values at the points in X
%--------------------------------------------------------------------------

    [mX,nX] = size(X); %dimensions of the points where function value should be predicted
    [mS,nS] = size(S); %dimesnions of sample sites taken so far 
    if nX ~= nS %check that both matrices are of the right shape
        X = X';
        [mX,nX] = size(X);
    end
    
    R = zeros(mX,mS); %compute pairwise distances of points in X and S
    for ii = 1 : mX
        for jj = 1 : mS
            R(ii,jj) = norm(X(ii,:)-S(jj,:));
        end
    end
    
    if strcmp(flag,'CB') %cubic RBF
        Phi = R.^3;
    elseif strcmp(flag,'TPS') %thin plate spline RBF
        Phi = (R.^2).*log(R);
        Phi(isnan(Phi)) = 0;
    elseif strcmp(flag,'LN') %linear RBF
        Phi = R;
    elseif strcmp(flag,'GA') %gaussian RBF
        Phi = exp(-(R.^2)./(rho^2));
    elseif strcmp(flag,'MQ') % multiquadric
        Phi = sqrt(((R./rho).^2)+1); 
    elseif strcmp(flag,'IMQ') % invers multiquadric
        Phi = 1./sqrt(((R./rho).^2)+1); 
    end
        
    Yest1 = Phi*lambda; %first part of response surface
    Yest2 = [X,ones(mX,1)]*gamma; %optional polynomial tail
    Yest  = Yest1+Yest2; %predicted function value
end
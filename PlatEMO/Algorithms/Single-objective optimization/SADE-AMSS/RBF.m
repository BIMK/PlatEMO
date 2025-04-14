function [lambda, gamma]=RBF(S,Y,flag)
%function [lambda, gamma, rho]=RBF(S,Y,flag)   % For Gaussian basis function

%computes the parameters of the radial basis function interpolant
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
%Input:
%S - Sample site matrix (m x d), m=number of sample sites, d=dimension
%Y - Objective function values corresponding to points in S
%flag - string determining what type of RBF model should be used
%
%Output:
%lambda, gamma - vectors with RBF model parameters
%--------------------------------------------------------------------------


[m,n]=size(S);
P=[S,ones(m,1)];

R=zeros(m,m);
R = pdist2(S,S,'euclidean');
%rho=2;
%Sd=real(sqrt(S.^2*ones(size(S'))+ones(size(S))*(S').^2-2*S*(S')));    % For Gaussian basis function
%rho=max(max(Sd))/(n*m)^(1/n);     % For Gaussian basis function

if strcmp(flag,'cubic') %cubic RBF
    Phi= R.^3;
elseif strcmp(flag,'TPS') %thin plate spline RBF
    R(R==0)=1;
    Phi=R.^2.*log(R);
elseif strcmp(flag,'linear') %linear RBF
    Phi=R;
elseif strcmp(flag,'Gaussian')
    Phi=exp(-R.^2/rho^2);
elseif strcmp(flag,'multiquad')
    Phi=sqrt((R.^2 + rho^2));
elseif strcmp(flag,'invmultiquad')
    Phi=1./sqrt(R.^2 + rho^2);
end

A=[Phi,P;P', zeros(n+1,n+1)];
RHS=[Y;zeros(n+1,1)];
params=A\RHS;
lambda=params(1:m);
gamma=params(m+1:end);

end %function
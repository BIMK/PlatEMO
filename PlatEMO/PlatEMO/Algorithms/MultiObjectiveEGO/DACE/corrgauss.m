function  [r, dr] = corrgauss(theta, d)
%CORRGAUSS  Gaussian correlation function,
%
%           n
%   r_i = prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
%          j=1
%
% If length(theta) = 1, then the model is isotropic:
% all  theta_j = theta .
%
% Call:    r = corrgauss(theta, d)
%          [r, dr] = corrgauss(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 
% hbn@imm.dtu.dk  
% Last update June 2, 2002
%
% make the code faster for MATLAB2016b
% zhandawei@hust.edu.cn

[m, n] = size(d);  % number of differences and dimension of data
td = d.^2 .* (repmat(-theta(:)',m,1));
r = exp(sum(td, 2));

if  nargout > 1
  dr = (ones(m,1)*(-2*theta(:).')) .* d .* (r*ones(1,n));
end
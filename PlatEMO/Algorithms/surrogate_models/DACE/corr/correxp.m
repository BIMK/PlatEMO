function  [r, dr] = correxp(theta, d)
%CORREXP  Exponential correlation function
%
%           n
%   r_i = prod exp(-theta_j * |d_ij|)
%          j=1
%
% If length(theta) = 1, then the model is isotropic: 
% theta_j = theta(1), j=1,...,n
%
% Call:    r = correxp(theta, d)
%          [r, dr] = correxp(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% hbn@imm.dtu.dk  
% Last update April 12, 2002

[m n] = size(d);  % number of differences and dimension of data
lt = length(theta);
if  lt == 1,  theta = repmat(theta,1,n);
elseif  lt ~= n
  error(sprintf('Length of theta must be 1 or %d',n))
else
  theta = theta(:).';
end

td = abs(d) .* repmat(-theta, m, 1);
r = exp(sum(td,2));

if  nargout > 1
  dr = repmat(-theta,m,1) .* sign(d) .* repmat(r,1,n);
end
function  [r, dr] = correxpg(theta, d)
%CORREXPG  General exponential correlation function
%
%           n
%   r_i = prod exp(-theta_j * d_ij^theta_n+1)
%          j=1
%
% If n > 1 and length(theta) = 2, then the model is isotropic: 
% theta_j = theta(1), j=1,...,n;  theta_(n+1) = theta(2) 
%
% Call:    r = correxpg(theta, d)
%          [r, dr] = correxpg(theta, d)
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
if  n > 1 & lt == 2
  theta = [repmat(theta(1),1,n)  theta(2)];
elseif  lt ~= n+1
  error(sprintf('Length of theta must be 2 or %d',n+1))
else
  theta = theta(:).';
end

pow = theta(end);   tt = repmat(-theta(1:n), m, 1);
td = abs(d).^pow .* tt;
r = exp(sum(td,2));

if  nargout > 1
  dr = pow  * tt .* sign(d) .* (abs(d) .^ (pow-1)) .* repmat(r,1,n);
end

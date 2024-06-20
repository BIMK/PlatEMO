function  [r, dr] = corrspherical(theta, d)
%CORRSPHERICAL  Spherical correlation function,
%
%           n
%   r_i = prod max(0, 1 - 1.5(theta_j*d_ij) + .5(theta_j*d_ij)^3) ,  i = 1,...,m
%          j=1
%
% If length(theta) = 1, then the model is isotropic:
% all  theta_j = theta .
%
% Call:    r = corrspherical(theta, d)
%          [r, dr] = corrspherical(theta, d)
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
if  length(theta) == 1
  theta = repmat(theta,1,n);
elseif  length(theta) ~= n
  error(sprintf('Length of theta must be 1 or %d',n))
else
  theta = theta(:).';
end
td = min(abs(d) .* repmat(theta,m,1), 1);
ss = 1 - td .* (1.5 - .5*td.^2);
r = prod(ss, 2);

if  nargout > 1
  dr = zeros(m,n);
  for  j = 1 : n
    dd = 1.5*theta(j) * sign(d(:,j)).*(td(:,j).^2 - 1);
    dr(:,j) = prod(ss(:,[1:j-1 j+1:n]),2) .* dd;
  end
end  
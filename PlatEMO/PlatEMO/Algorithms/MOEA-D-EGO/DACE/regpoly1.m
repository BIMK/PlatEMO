function  [f, df] = regpoly1(S)
%REGPOLY1  First order polynomial regression function
%
% Call:    f = regpoly1(S)
%          [f, df] = regpoly1(S)
%
% S : m*n matrix with design sites
% f = [1  s]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update April 12, 2002

[m n] = size(S);
f = [ones(m,1)  S];
if  nargout > 1
  df = [zeros(n,1) eye(n)];
end
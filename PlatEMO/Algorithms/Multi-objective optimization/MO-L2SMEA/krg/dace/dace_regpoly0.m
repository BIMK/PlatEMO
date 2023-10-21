function  [f, df] = dace_regpoly0(S)
%Zero order polynomial regression function
%
% Call:    f = dace_regpoly0(S)
%          [f, df] = dace_regpoly0(S)
%
% S  : m*n matrix with design sites
% f  : ones(m,1)
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update  April 12, 2002

[m n] = size(S);
f = ones(m,1);
if  nargout > 1
  df = zeros(n,1);
end

return

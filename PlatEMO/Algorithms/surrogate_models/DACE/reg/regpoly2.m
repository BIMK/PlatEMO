function  [f, df] = regpoly2(S)
%REGPOLY2  Second order polynomial regression function
% Call:    f = regpoly2(S)
%          [f, df] = regpoly2(S)
%
% S : m*n matrix with design sites
% f =  [1 S S(:,1)*S S(:,2)S(:,2:n) ... S(:,n)^2]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update September 4, 2002

[m n] = size(S);
nn = (n+1)*(n+2)/2;  % Number of columns in f  
% Compute  f
f = [ones(m,1) S zeros(m,nn-n-1)];
j = n+1;   q = n;
for  k = 1 : n
  f(:,j+(1:q)) = repmat(S(:,k),1,q) .* S(:,k:n);
  j = j+q;   q = q-1;
end

if  nargout > 1
  df = [zeros(n,1)  eye(n)  zeros(n,nn-n-1)];
  j = n+1;   q = n; 
  for  k = 1 : n
    df(k,j+(1:q)) = [2*S(1,k) S(1,k+1:n)];
    for i = 1 : n-k,  df(k+i,j+1+i) = S(1,k); end
    j = j+q;   q = q-1;
  end
end 
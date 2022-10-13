% KERNEL - Kernel function evaluation  
%
% Syntax: K = kernel(X,Y,type,scale)
%
%      X: (N,Mx) dimensional matrix
%      Y: (N,My) dimensional matrix
%      K: (Mx,My) dimensional matrix 
%   type: kernel type
%           1: linear kernel        X'*Y
%         2-4: polynomial kernel    (scale*X'*Y + 1)^type
%           5: Gaussian kernel with variance 1/(2*scale)
%  scale: kernel scale
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function K = kernel(X,Y,type,scale)

global kernel_evals;			% kernel evaluations

K = X'*Y;
[N,Mx] = size(X);
[N,My] = size(Y);
if ((type > 1) & (type < 5))
   K = (K*scale+1).^type;
elseif (type == 5)
   K = 2*K;
   K = K - sum(X.^2,1)'*ones(1,My);
   K = K - ones(Mx,1)*sum(Y.^2,1);
   K = exp(K/(2*scale));
end;
kernel_evals = kernel_evals + Mx*My;


% SVMEVAL - Evaluates a support vector machine at the given data points.
%
% Syntax: [f,K] = svmeval(X,a,b,ind,X_mer,y_mer,type,scale)
%         (evaluates the given SVM at the data points contained in X)
%
%         [f,K] = svmeval(X)
%         (evaluates the SVM in memory at the data points contained in X)
%
%      f: SVM output for the evaluation vectors
%      K: kernel matrix containing dot products in feature space between
%         the margin and error vectors (rows of K) and the column vectors in X
%         (columns of K)
%      X: matrix of evaluation vectors stored columnwise
%      a: alpha coefficients
%      b: bias
%    ind: cell array containing indices of margin, error and reserve vectors
%         ind{1}: indices of margin vectors
%         ind{2}: indices of error vectors
%         ind{3}: indices of reserve vectors
%  X_mer: matrix of margin, error and reserve vectors stored columnwise
%  y_mer: column vector of class labels (-1/+1) for margin, error and reserve vectors
%   type: kernel type
%           1: linear kernel        K(x,y) = x'*y
%         2-4: polynomial kernel    K(x,y) = (scale*x'*y + 1)^type
%           5: Gaussian kernel with variance 1/(2*scale)
%  scale: kernel scale
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [f,K] = svmeval(X_eval, varargin)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

if (nargin == 8)
   
   % define arguments
   a     = varargin{1};
   b     = varargin{2};
   ind   = varargin{3};
   X     = varargin{4};
   y     = varargin{5};
   type  = varargin{6};
   scale = varargin{7};
   
else
   
   % define global variables
   global a;
   global b;
   global ind;
   global scale;
   global type;
   global X;           
   global y;          
   
end;

% evaluate the SVM

% find all of the nonzero coefficients
% (note: when performing kernel perturbation, ind{MARGIN} and ind{ERROR}
%  do not necessarily identify all of the nonzero coefficients)
indu = find(a(ind{UNLEARNED}) > 0);
indu = ind{UNLEARNED}(indu);
indr = find(a(ind{RESERVE}) > 0);
indr = ind{RESERVE}(indr);
indme = [ind{MARGIN} ind{ERROR}];

K = [];
f = b;
if (length(indme) > 0)
   K = kernel(X(:,indme),X_eval,type,scale);
   f = f + K'*(y(indme).*a(indme));
end;
if (length(indu) > 0)
   f = f + kernel(X(:,indu),X_eval,type,scale)'*(y(indu).*a(indu));
end;
if (length(indr) > 0)
   f = f + kernel(X(:,indr),X_eval,type,scale)'*(y(indr).*a(indr));
end;


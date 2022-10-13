% SVMLOO - Computes the exact leave-one-out estimate of the error rate
%          for a two-class SVM using Cauwenbergh's method.
%
% Syntax: [loo_est,conf_matrix,g_loo] = svmloo;
%         (evaluates the current SVM in memory)
%       
%         [loo_est,conf_matrix,g_loo] = svmloo(X,y,a,b,g,ind,type,scale,Rs,Q)
%         (evaluates the given SVM)
%
%     loo_est: number of leave-one-out error errors
% conf_matrix: confusion matrix
%       g_loo: the resulting g for each example after the example is unlearned
%           X: matrix of training vectors stored columnwise
%           y: column vector of class labels (-1/+1) for training vectors
%           a: alpha coefficients
%           b: bias
%           g: partial derivatives of cost function w.r.t. alpha coefficients
%         ind: cell array containing indices of margin, error and reserve vectors
%         		ind{1}: indices of margin vectors
%         		ind{2}: indices of error vectors
%         		ind{3}: indices of reserve vectors
%        type: kernel type
%                1: linear kernel        X'*Y
%              2-4: polynomial kernel    (scale*X'*Y + 1)^type
%                5: Gaussian kernel with variance 1/(2*scale)
%       scale: kernel scale
%          Rs: inverse of extended kernel matrix for margin vectors
%           Q: extended kernel matrix for all vectors
%      g_flag: flag indicating whether or not to compute g_loo for the error vectors with g < -1
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [loo_est,conf_matrix,g_loo] = svmloo(varargin)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

if (nargin == 0)

   % define global variables 
   global a; 
   global b; 
   global g;                           
   global ind;
   global Q;
   global Rs;   
   global scale;
   global type;
   global X;
   global y;
   
else   
   
   % define arguments
   X = varargin{1};
   y = varargin{2};
   a = varargin{3};
   b = varargin{4};
   g = varargin{5};
   ind = varargin{6};
   type = varargin{7};
   scale = varargin{8};
   Rs = varargin{9};
   Q = varargin{10};
   g_flag = varargin(11);
   
end;

% if the user wants g_loo, make sure to compute g_loo for the error vectors with initial g < -1.
% if we only care about the error rate, we don't need to unlearn these examples because they are
% guaranteed to be classified incorrectly.
if (nargout == 3)
   g_flag = 1;
else
   g_flag = 0;
end;

% initialize variables
num_MVs = length(ind{MARGIN});      % number of margin vectors
loo_est = 0;                        % number of leave-one-out errors
a_orig = a;                         % original value of a
b_orig = b;                         % original value of b
Rs_orig = Rs;                       % original value of Rs
Q_orig = Q;                         % original value of Q                                  
g_orig = g;                         % original value of g
num_MVs_orig = num_MVs;             % original value of num_MVs
ind_orig = ind;                     % original value of ind
g_loo = g;

% initialize confusion matrix
conf_matrix = zeros(2,3);
if (length(ind{RESERVE}) > 0)
	conf_matrix(1,1) = sum(y(ind{RESERVE}) == 1);
	conf_matrix(2,2) = sum(y(ind{RESERVE}) == -1);
end;

% begin leave-one-out estimation
ind_loo = [ind{MARGIN} ind{ERROR}];
num_tested = 1;
disp('Beginning LOO error rate estimation.');
for i = 1:length(ind_loo)
   
   % select example to unlearn
   indc = ind_loo(i);
   
   % unlearn example
   if ((g(indc) >= -1) | (g_flag))
      unlearn(indc);
      g_loo(indc) = g(indc);   
   end;
   
   if (mod(num_tested,50) == 0)
      s = sprintf('Unlearned and tested %d examples.',num_tested);
      disp(s);
   end;
   num_tested = num_tested + 1;
   
   % check to see if the example is now misclassified and record results
   loo_est = loo_est + (g(indc) < -1);
   if (y(indc) == 1)
      j = 1 + (g(indc) <= -1) + (g(indc) == -1);
      conf_matrix(1,j) = conf_matrix(1,j) + 1;
   else
      j = 1 + (g(indc) >= -1) + (g(indc) == -1);
      conf_matrix(2,j) = conf_matrix(2,j) + 1;
   end;
      
   % reset to original state prior to unlearning example
   a = a_orig;
   b = b_orig;
   Rs = Rs_orig;
   Q = Q_orig;                                 
   g = g_orig;
   ind = ind_orig;
   num_MVs = num_MVs_orig;
   
end;
if (mod(num_tested-1,50) ~= 0)
   s = sprintf('Unlearned and tested %d examples.',num_tested-1);
   disp(s);
end;
s = sprintf('Process complete!\n');
disp(s);

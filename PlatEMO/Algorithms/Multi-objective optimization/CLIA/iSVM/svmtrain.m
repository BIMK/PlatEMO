% SVMTRAIN - Trains a support vector machine incrementally
%            using the L1 soft margin approach developed by
%            Cauwenberghs for two-class problems.
%
% Syntax: [a,b,g,ind,uind,X_mer,y_mer,Rs,Q] = svmtrain(X,y,C,type,scale)
%         [a,b,g,ind,uind,X_mer,y_mer,Rs,Q] = svmtrain(X,y,C,type,scale,uind)
%         (trains a new SVM on the given examples)
%
%         [a,b,g,ind,uind,X_mer,y_mer,Rs,Q] = svmtrain(X,y,C)
%         [a,b,g,ind,uind,X_mer,y_mar,Rs,Q] = svmtrain(X,y,C,uind)
%         (trains the current SVM in memory on the given examples)
%
%      a: alpha coefficients
%      b: bias
%      g: partial derivatives of cost function w.r.t. alpha coefficients
%    ind: cell array containing indices of margin, error and reserve vectors
%         ind{1}: indices of margin vectors
%         ind{2}: indices of error vectors
%         ind{3}: indices of reserve vectors
%   uind: column vector of user-defined example indices (used for unlearning specified examples)
%  X_mer: matrix of margin, error and reserve vectors stored columnwise
%  y_mer: column vector of class labels (-1/+1) for margin, error and reserve vectors
%     Rs: inverse of extended kernel matrix for margin vectors
%      Q: extended kernel matrix for all vectors
%      X: matrix of training vectors stored columnwise
%      y: column vector of class labels (-1/+1) for training vectors
%      C: soft-margin regularization parameter(s)
%         dimensionality of C       assumption
%         1-dimensional vector      universal regularization parameter
%         2-dimensional vector      class-conditional regularization parameters (-1/+1)
%         n-dimensional vector      regularization parameter per example
%         (where n = # of examples)
%   type: kernel type
%           1: linear kernel        K(x,y) = x'*y
%         2-4: polynomial kernel    K(x,y) = (scale*x'*y + 1)^type
%           5: Gaussian kernel with variance 1/(2*scale)
%  scale: kernel scale
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [a,b,g,ind,uind,X,y,Rs,Q] = svmtrain(X_new,y_new,C_new,varargin)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

% create a vector containing the regularization parameter 
% for each example if necessary
if (length(C_new) == 1)              % same regularization parameter for all examples
   C_new = C_new*ones(size(y_new));    
elseif (length(C_new) == 2)          % class-conditional regularization parameters
   flags = (y_new == -1);
   C_new = C_new(1)*flags + C_new(2)*(~flags);
end;

if (nargin >= 5)
   
   % define arguments      
   type_new = varargin{1};
   scale_new = varargin{2};
   if (nargin == 6)
      uind_new = varargin(3);
   else
      uind_new = zeros(size(y_new));
   end;
         
   new_model = 1;   
else
   
   % define arguments
   if (nargin == 4)
      uind_new = varargin(1);
   else
      uind_new = zeros(size(y_new));
   end;
      
   new_model = 0;   
end;

% define global variables 
global a;                     % alpha coefficients
global b;                     % bias
global C;                     % regularization parameters 
global deps;                  % jitter factor in kernel matrix
global g;                     % partial derivatives of cost function w.r.t. alpha coefficients
global ind;                   % cell array containing indices of margin, error, reserve and unlearned vectors
global kernel_evals;          % kernel evaluations
global max_reserve_vectors;   % maximum number of reserve vectors stored
global perturbations;         % number of perturbations
global Q;                     % extended kernel matrix for all vectors
global Rs;                    % inverse of extended kernel matrix for margin vectors   
global scale;                 % kernel scale
global type;                  % kernel type
global uind;                  % user-defined example indices
global X;                     % matrix of margin, error, reserve and unlearned vectors stored columnwise
global y;                     % column vector of class labels (-1/+1) for margin, error, reserve and unlearned vectors

% initialize variables
deps = 1e-3;
max_reserve_vectors = 3000;    

if (new_model)

   num_examples = size(X_new,2);       
   
   a = zeros(num_examples,1);          
   b = 0;                              
   C = C_new;                          
   g = zeros(num_examples,1);
   ind = cell(4,1);
   ind{UNLEARNED} = 1:num_examples;
   kernel_evals = 0;
   perturbations = 0;
   Q = y_new';
   Rs = Inf;
   scale = scale_new;
   type = type_new;
   uind = uind_new;
   X = X_new;                          
   y = y_new;

else
   
   num_examples = size(X,2);
   num_new_examples = size(X_new,2);
   
   a = [a ; zeros(num_new_examples,1)];
   C = [C ; C_new];
   g = [g ; zeros(num_new_examples,1)];
   ind{UNLEARNED} = (1:num_new_examples) + num_examples;
   
   % assumes currently that there are no duplicate examples in the data - may not necessarily be true!
   Q_new = [y_new' ; (y(ind{MARGIN})*y_new').*kernel(X(:,ind{MARGIN}),X_new,type,scale)];
   
   Q = [Q Q_new];  
   uind = [uind ; uind_new];
   X = [X X_new];
   y = [y ; y_new];
   
   num_examples = num_examples + num_new_examples;
   
end;   
   
% begin incremental learning - enforce all constraints on each iteration
num_learned = 1;
disp('Beginning training.');
while (any(ind{UNLEARNED}))
   
   % randomly select example
   i = round(rand*(length(ind{UNLEARNED})-1)) + 1;
   indc = ind{UNLEARNED}(i);
%  indc = ind{UNLEARNED}(1);

	% learn example
   learn(indc,1);
   
   if (mod(num_learned,50) == 0)
      s = sprintf('Learned %d examples.',num_learned);
      disp(s);
   end;
   num_learned = num_learned + 1;
   
end;
if (mod(num_learned-1,50) ~= 0)
   s = sprintf('Learned %d examples.',num_learned-1);
   disp(s);
end;
disp('Training complete!');

% begin incremental learning - perform multiple passes through the data 
% until all of the examples are learned
%while (any(ind{UNLEARNED}))
%   while (any(ind{UNLEARNED}))
%   
%      % select example
%      indc = ind{UNLEARNED}(1);
%      
%      % learn example
%      s = sprintf('\nLearning example %d...',indc);
%      disp(s);
%      learn(indc,0);
%   
%   end;
%
%   % check to see if any reserve vectors are incorrectly classified
%   % if so, change their status to unlearned
%   ind_temp = find(g(ind{RESERVE}) < 0);
%   [ind{RESERVE},ind{UNLEARNED}] = move_ind(ind{RESERVE},ind{UNLEARNED},ind{RESERVE}(ind_temp));
%   
%end;

% remove all but the closest reserve vectors from the dataset if necessary
if (length(ind{RESERVE}) == max_reserve_vectors)
   ind_keep = [ind{MARGIN} ind{ERROR} ind{RESERVE}];
   a = a(ind_keep);
   g = g(ind_keep);
   Q = Q(:,ind_keep);   
   uind = uind(ind_keep);
   X = X(:,ind_keep);
   y = y(ind_keep);
   ind{MARGIN} = 1:length(ind{MARGIN});
   ind{ERROR} = length(ind{MARGIN}) + (1:length(ind{ERROR}));
   ind{RESERVE} = length(ind{MARGIN}) + length(ind{ERROR}) + (1:length(ind{RESERVE}));
end;

% summary statistics
s = sprintf('\nMargin vectors:\t\t%d',length(ind{MARGIN}));
disp(s);
s = sprintf('Error vectors:\t\t%d',length(ind{ERROR}));
disp(s);
s = sprintf('Reserve vectors:\t%d',length(ind{RESERVE}));
disp(s);
s = sprintf('Kernel evaluations:\t%d\n',kevals);
disp(s);

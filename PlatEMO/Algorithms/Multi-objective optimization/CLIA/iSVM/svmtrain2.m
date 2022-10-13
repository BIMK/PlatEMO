
% SVMTRAIN2 - Trains a support vector machine in batch mode
%             using the L1 soft margin approach developed by
%             Diehl and Cauwenberghs for two-class problems.
%
% Syntax: [a,b,g,ind,uind,X_mer,y_mer,Rs,Q] = svmtrain2(X,y,C,type,scale)
%         [a,b,g,ind,uind,X_mer,y_mer,Rs,Q] = svmtrain2(X,y,C,type,scale,uind)
%         (trains a new SVM on the given examples)
%
%         [a,b,g,ind,uind,X_mer,y_mer,Rs,Q] = svmtrain2(X,y,C)
%         [a,b,g,ind,uind,X_mer,y_mer,Rs,Q] = svmtrain2(X,y,C,uind)
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

function [a,b,g,ind,uind,X,y,Rs,Q] = svmtrain2(X_new,y_new,C_new,varargin)

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
deps = 1e-5;
max_reserve_vectors = 3000;    

if (new_model)

   num_examples = size(X_new,2);       
   
   a = zeros(num_examples,1);          
   b = 0;                              
   C = C_new;                          
   g = -ones(num_examples,1);
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
   
   g_new = y_new.*svmeval(X_new) - 1;
   g = [g ; g_new];
   flag = (g_new > 0);
   indr = find(flag)';
   indu = find(~flag)';
   move_indr([],indr+num_examples);
   ind{UNLEARNED} = indu+num_examples;

   % assumes currently that there are no duplicate examples in the data - may not necessarily be true!
   if (length(ind{MARGIN}) > 0)
      Q_new = [y_new' ; (y(ind{MARGIN})*y_new').*kernel(X(:,ind{MARGIN}),X_new,type,scale)];
   else
      Q_new = y_new';
   end;
   Q = [Q Q_new];   
   uind = [uind ; uind_new];
   X = [X X_new];
   y = [y ; y_new];
   
   num_examples = num_examples + num_new_examples;
   
end;   
   
% perturbation coefficients
lambda = zeros(num_examples,1);
lambda(ind{UNLEARNED}) = C(ind{UNLEARNED});

% cached results
disp('Precomputing results for unlearned vectors.');
SQl = ((y*y(ind{UNLEARNED})').*kernel(X,X(:,ind{UNLEARNED}),type,scale))*lambda(ind{UNLEARNED}) + deps*lambda;
Syl = y(ind{UNLEARNED})'*lambda(ind{UNLEARNED});
    
p = 0;
num_MVs = length(ind{MARGIN});
num_learned = 0;
disp('Beginning training.');
i = 0;
while ((p < 1) | (length(ind{UNLEARNED}) > 0))
   
   if (p < 1)
   
      perturbations = perturbations + 1;
   
      % compute beta and gamma
      if (num_MVs > 0)
      
         v = zeros(num_MVs+1,1);
         if (p < 1-eps)
            v(1) = -Syl - sum(y.*a)/(1-p);
         else
            v(1) = -Syl;
         end;
         v(2:num_MVs+1) = -SQl(ind{MARGIN});
         beta = Rs*v;
         gamma = zeros(size(Q,2),1);
         ind_temp = [ind{ERROR} ind{RESERVE} ind{UNLEARNED}];
         if (length(ind_temp) > 0)
            gamma(ind_temp) = Q(:,ind_temp)'*beta + SQl(ind_temp);
         end;
         
      else
      
         beta = 0;
         gamma = SQl;
                  
      end;
        
      % minimum acceptable parameter change
      [min_delta_p,indss,cstatus,nstatus] = min_delta_p_s(p,gamma,beta,lambda);
   
      % update a, b, g and p
      if (length(ind{UNLEARNED}) > 0)
         a(ind{UNLEARNED}) = a(ind{UNLEARNED}) + lambda(ind{UNLEARNED})*min_delta_p;
      end;
      if (num_MVs > 0)
         a(ind{MARGIN}) = a(ind{MARGIN}) + beta(2:num_MVs+1)*min_delta_p;
      end;   
      b = b + beta(1)*min_delta_p;
      g = g + gamma*min_delta_p;
      p = p + min_delta_p;
   
      % perform bookkeeping         
      indco = bookkeeping(indss,cstatus,nstatus);
   
      % update Rs and Q if necessary
      if (nstatus == MARGIN)
              
         num_MVs = num_MVs + 1;
         if (num_MVs > 1)
               
            % compute beta and gamma for indss            
            beta = -Rs*Q(:,indss);
            gamma = kernel(X(:,indss),X(:,indss),type,scale) + deps + Q(:,indss)'*beta;
            
         end;
            
         % expand Rs and Q
         updateRQ(beta,gamma,indss);
      
      elseif (cstatus == MARGIN)      
              
         % compress Rs and Q      
         num_MVs = num_MVs - 1;
         updateRQ(indco);
         
      end;         
   
      % update SQl and Syl when the current status of indss is UNLEARNED
      if (cstatus == UNLEARNED)
         num_learned = num_learned + 1;
         if (nstatus == MARGIN)
            SQl = SQl - Q(num_MVs+1,:)'*lambda(indss);  
         else
            SQl = SQl - ((y*y(indss)).*kernel(X,X(:,indss),type,scale))*lambda(indss);
            SQl(indss) = SQl(indss) - deps*lambda(indss);
         end;
         Syl = Syl - y(indss)*lambda(indss);
         
         if (mod(num_learned,50) == 0)
            s = sprintf('Learned %d examples.',num_learned);
            disp(s);
         end;
      end;
   
   else
      
      % label the remaining unlearned vectors as error vectors
      num_learned = num_learned + length(ind{UNLEARNED});
      a(ind{UNLEARNED}) = C(ind{UNLEARNED});
      [ind{UNLEARNED},ind{ERROR}] = move_ind(ind{UNLEARNED},ind{ERROR},ind{UNLEARNED});
      
      if (mod(num_learned,50) == 0)
         s = sprintf('Learned %d examples.',num_learned);
         disp(s);
      end;
      
   end;
   
end;
if (mod(num_learned,50) ~= 0)
   s = sprintf('Learned %d examples.',num_learned);
   disp(s);
end;
disp('Training complete!');

% set g(ind{MARGIN}) to zero
g(ind{MARGIN}) = 0;

% remove all but the closest reserve vectors from the dataset if necessary
if (length(ind{RESERVE}) == max_reserve_vectors)
   ind_keep = [ind{MARGIN} ind{ERROR} ind{RESERVE}];
   X = X(:,ind_keep);
   y = y(ind_keep);
   a = a(ind_keep);
   g = g(ind_keep);
   Q = Q(:,ind_keep);   
   uind = uind(ind_keep);
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

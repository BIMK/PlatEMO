% PERTURBC - Perturbs the current solution to the solution valid for the
%            given regularization parameters.  
%
% Syntax: [a,b,g,ind,X_mer,y_mer,Rs,Q] = perturbc(C)
%
%      a: alpha coefficients
%      b: bias
%      g: partial derivatives of cost function w.r.t. alpha coefficients
%    ind: cell array containing indices of margin, error and reserve vectors
%         ind{1}: indices of margin vectors
%         ind{2}: indices of error vectors
%         ind{3}: indices of reserve vectors
%  X_mer: matrix of margin, error and reserve vectors stored columnwise
%  y_mer: column vector of class labels (-1/+1) for margin, error and reserve vectors
%     Rs: inverse of extended kernel matrix for margin vectors
%      Q: extended kernel matrix for all vectors
%      C: soft-margin regularization parameter(s)
%         dimensionality of C       assumption
%         1-dimensional vector      universal regularization parameter
%         2-dimensional vector      class-conditional regularization parameters (-1/+1)
%         n-dimensional vector      regularization parameter per example
%         (where n = # of examples)
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [a,b,g,ind,X,y,Rs,Q] = perturbc(C_new)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

% define global variables
global a;               % alpha coefficients
global b;               % bias
global C;               % regularization parameters 
global deps;            % jitter factor in kernel matrix
global g;               % partial derivatives of cost function w.r.t. alpha coefficients
global ind;             % cell array containing indices of margin, error, reserve and unlearned vectors
global perturbations;   % number of perturbations
global Q;               % extended kernel matrix for all vectors
global Rs;              % inverse of extended kernel matrix for margin vectors   
global scale;           % kernel scale
global type;            % kernel type
global X;               % matrix of margin, error, reserve and unlearned vectors stored columnwise
global y;               % column vector of class labels (-1/+1) for margin, error, reserve and unlearned vectors

kernel_evals_begin = kevals;

% create a vector containing the regularization parameter 
% for each example if necessary
if (length(C_new) == 1)              % same regularization parameter for all examples
   C_new = C_new*ones(size(y));    
elseif (length(C_new) == 2)          % class-conditional regularization parameters
   flags = (y == -1);
   C_new = C_new(1)*flags + C_new(2)*(~flags);
end;

% compute the regularization sensitivities
lambda = C_new-C;

% if there are no error vectors initially...
if (length(ind{ERROR}) == 0)
   
   % find all the examples that have changing regularization parameters   
   inde = find(lambda ~= 0);
   
   % find the subset of the above examples that could become error vectors
   delta_p = (a(inde)-C(inde))./lambda(inde);
   i = find(delta_p > 0);
   
   % determine the minimum acceptable change in p and adjust the regularization parameters
   p = min([delta_p(i) ; 1]);
   C = C + lambda*p;
   
   % if one example becomes an error vector, perform the necessary bookkeeping
   if (p < 1)
      i = find(delta_p == p);
      indco = bookkeeping(inde(i),MARGIN,ERROR);
      updateRQ(indco);
   end;
   
else
   p = 0;
end;

% if there are error vectors to adjust...
if (p < 1)
   
   % compute sum{k in E} Qik lambda k and sum{k in E} yk lambda k
   SQl = ((y*y(ind{ERROR})').*kernel(X,X(:,ind{ERROR}),type,scale))*lambda(ind{ERROR});
   SQl(ind{ERROR}) = SQl(ind{ERROR}) + deps*lambda(ind{ERROR});
   Syl = y(ind{ERROR})'*lambda(ind{ERROR});
   
end;   
      
s = sprintf('p = %.2f',p);
disp(s);

% change the regularization parameters incrementally
disp_p_delta = 0.2;
disp_p_count = 1;
num_MVs = length(ind{MARGIN});
perturbations = 0;
while (p < 1) 
   
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
   [min_delta_p,indss,cstatus,nstatus] = min_delta_p_c(p,gamma,beta,lambda);
   
   % update a, b, g and p
   if (length(ind{ERROR}) > 0)
      a(ind{ERROR}) = a(ind{ERROR}) + lambda(ind{ERROR})*min_delta_p;
   end;
   if (num_MVs > 0)
      a(ind{MARGIN}) = a(ind{MARGIN}) + beta(2:num_MVs+1)*min_delta_p;
   end;   
   b = b + beta(1)*min_delta_p;
   g = g + gamma*min_delta_p;
   p = p + min_delta_p;
   C = C + lambda*min_delta_p;
   
   % perform bookkeeping         
   indco = bookkeeping(indss,cstatus,nstatus);
   
   % update SQl and Syl when the status of indss changes from MARGIN to ERROR
   if ((cstatus == MARGIN) & (nstatus == ERROR))
      SQl = SQl + Q(indco,:)'*lambda(indss);
      Syl = Syl + y(indss)*lambda(indss);      
   end;
   
   % set g(ind{MARGIN}) to zero
   g(ind{MARGIN}) = 0;
   
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
   
   % update SQl and Syl when the status of indss changes from ERROR to MARGIN
   if ((cstatus == ERROR) & (nstatus == MARGIN))
      SQl = SQl - Q(num_MVs+1,:)'*lambda(indss);
      Syl = Syl - y(indss)*lambda(indss);
   end;
   
   if (p >= disp_p_delta*disp_p_count)
      disp_p_count = disp_p_count + 1;
      s = sprintf('p = %.2f',p);
      disp(s);
   end;
  
end;
disp('Perturbation complete!');

% summary statistics
s = sprintf('\nMargin vectors:\t\t%d',length(ind{MARGIN}));
disp(s);
s = sprintf('Error vectors:\t\t%d',length(ind{ERROR}));
disp(s);
s = sprintf('Reserve vectors:\t%d',length(ind{RESERVE}));
disp(s);
s = sprintf('Kernel evaluations:\t%d\n',-kernel_evals_begin+kevals);
disp(s);

   

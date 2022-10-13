% PERTURBK - Perturbs the current solution to the solution valid for the
%            given kernel parameter.
%
% Syntax: [a,b,g,ind,X_mer,y_mer,Rs,Q] = perturbk(scale)
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
%  scale: kernel scale
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [a,b,g,ind,X,y,Rs,Q] = perturbk(new_scale)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

% define global variables
global a;             % alpha coefficients
global b;             % bias
global C;             % regularization parameters 
global deps;          % jitter factor in kernel matrix
global g;             % partial derivatives of cost function w.r.t. alpha coefficients
global ind;           % cell array containing indices of margin, error, reserve and unlearned vectors
global num_unlearned; % number of unlearned vectors initially
global perturbations; % number of perturbations
global Q;             % extended kernel matrix for all vectors
global Rs;            % inverse of extended kernel matrix for margin vectors   
global scale;			 % kernel scale
global type;          % kernel type
global X;             % matrix of margin, error, reserve and unlearned vectors stored columnwise
global y;             % column vector of class labels (-1/+1) for margin, error, reserve and unlearned vectors

kernel_evals_begin = kevals;
num_examples = size(X,2);

% sum{k in U} Qik lambda k and sum{k in U} yk lambda k
SQl = zeros(num_examples,1);
Syl = 0;

% adjust kernel scale
scale = new_scale;

% recompute g for the margin, error and reserve vectors
inda = [ind{MARGIN} ind{ERROR} ind{RESERVE}];
[f,K] = svmeval(X(:,inda));
g(inda) = y(inda).*f + a(inda)*deps - 1;

% identify the unlearned vectors and compute their coefficient sensitivities
lambda = zeros(num_examples,1);

% find all error vectors with g >= 0
flag = (g(ind{ERROR}) >= 0);
i = find(flag);

if (length(i) > 0)

   % relabel error vectors with g >= 0 as unlearned
   ind{UNLEARNED} = [ind{UNLEARNED} ind{ERROR}(i)];

   % coefficient sensitivities
   lambda(ind{ERROR}(i)) = -a(ind{ERROR}(i));

   % update sums
   SQl(inda) = SQl(inda) + ((y(ind{ERROR}(i))*y(inda)').*K(length(ind{MARGIN})+i,:))'*lambda(ind{ERROR}(i));
   Syl = Syl + y(ind{ERROR}(i))'*lambda(ind{ERROR}(i));
         
   % keep remaining error vectors labeled as such
   ind{ERROR}(i) = [];
   
end;
   
% find all reserve vectors with g <= 0
flag = (g(ind{RESERVE}) <= 0);
i = find(flag);

if (length(i) > 0)

   % relabel reserve vectors with g <= 0 as unlabeled
   ind{UNLEARNED} = [ind{UNLEARNED} ind{RESERVE}(i)];

   % coefficient sensitivities
   lambda(ind{RESERVE}(i)) = C(ind{RESERVE}(i));

   % update sums
   SQl(inda) = SQl(inda) + ((y(inda)*y(ind{RESERVE}(i))').*kernel(X(:,inda),X(:,ind{RESERVE}(i)),type,scale))*lambda(ind{RESERVE}(i));
   Syl = Syl + y(ind{RESERVE}(i))'*lambda(ind{RESERVE}(i));
      
   % keep remaining reserve vectors labeled as such
   ind{RESERVE}(i) = [];
   
end;   
   
% find all margin vectors with g > 0
flag = (g(ind{MARGIN}) > 0);
i = find(flag);

if (length(i) > 0)

   % coefficient sensitivities
   lambda(ind{MARGIN}(i)) = -a(ind{MARGIN}(i));

   % update sums
   SQl(inda) = SQl(inda) + ((y(ind{MARGIN}(i))*y(inda)').*K(i,:))'*lambda(ind{MARGIN}(i));
   Syl = Syl + y(ind{MARGIN}(i))'*lambda(ind{MARGIN}(i));
         
end;   
   
% find all margin vectors with g <= 0
flag = ~flag;
i = find(flag);

if (length(i) > 0)

   % coefficient sensitivities
   lambda(ind{MARGIN}(i)) = C(ind{MARGIN}(i))-a(ind{MARGIN}(i));

   % update sums
   SQl(inda) = SQl(inda) + ((y(ind{MARGIN}(i))*y(inda)').*K(i,:))'*lambda(ind{MARGIN}(i));
   Syl = Syl + y(ind{MARGIN}(i))'*lambda(ind{MARGIN}(i));
      
end;   
   
% relabel margin vectors as unlearned
ind{UNLEARNED} = [ind{UNLEARNED} ind{MARGIN}];
ind{MARGIN} = [];

% add jitter factor
SQl(ind{UNLEARNED}) = SQl(ind{UNLEARNED}) + deps*lambda(ind{UNLEARNED});

% number of unlearned vectors initially
num_unlearned = length(ind{UNLEARNED});
s = sprintf('Number of unlearned vectors: %d',num_unlearned);
disp(s);

% reset Q and Rs
Q = Q(1,:);
Rs = Inf;

p_s = 0;
num_MVs = length(ind{MARGIN});
num_learned = 0;
perturbations = 0;
while ((length(ind{UNLEARNED}) > 0) | ((p_s < 1) & (length(ind{UNLEARNED}) == 0)))
   
   perturbations = perturbations + 1;
   
   % compute beta and gamma
   if (num_MVs > 0)
      
      v = zeros(num_MVs+1,1);
      if (p_s < 1-eps)
         v(1) = -Syl - sum(y.*a)/(1-p_s);
      else
         v(1) = -Syl;
      end;
      v(2:num_MVs+1) = -SQl(ind{MARGIN});
      beta = Rs*v;
      gamma = zeros(size(Q,2),1);
      ind_temp = [ind{ERROR} ind{RESERVE} ind{UNLEARNED}];
      gamma(ind_temp) = Q(:,ind_temp)'*beta + SQl(ind_temp);
      
   else
      
      beta = 0;
      gamma = SQl;
                  
   end;
        
   % minimum acceptable parameter change
   [min_dps,indss,cstatus,nstatus] = min_delta_p_s(p_s,gamma,beta,lambda);
   
   % update a, b, g and p_s
   if (length(ind{UNLEARNED}) > 0)
      a(ind{UNLEARNED}) = a(ind{UNLEARNED}) + lambda(ind{UNLEARNED})*min_dps;
   end;
   if (num_MVs > 0)
      a(ind{MARGIN}) = a(ind{MARGIN}) + beta(2:num_MVs+1)*min_dps;
   end;   
   b = b + beta(1)*min_dps;
   g = g + gamma*min_dps;
   p_s = p_s + min_dps;
   
   % perform bookkeeping         
   indco = bookkeeping(indss,cstatus,nstatus);
   
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


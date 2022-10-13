% LEARN - Increments the specified example into the current SVM solution.  
%         Assumes alpha_c = 0 initially.
%
% Syntax: nstatus = learn(indc,rflag)
%
% nstatus: new status for indc
%    indc: index of the example to learn
%   rflag: flag indicating whether or not to check if any reserve vectors
%          become margin vectors during learning
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function nstatus = learn(indc,rflag)

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
global ind;           % structure containing indices of margin, error, reserve and unlearned vectors
global perturbations; % number of perturbations
global Q;             % extended kernel matrix for all vectors
global Rs;            % inverse of extended kernel matrix for margin vectors   
global scale;         % kernel scale
global type;          % kernel type
global X;             % matrix of margin, error, reserve and unlearned vectors stored columnwise
global y;             % column vector of class labels (-1/+1) for margin, error, reserve and unlearned vectors

% compute g(indc) 
[f_c,K] = svmeval(X(:,indc));
g(indc) = y(indc)*f_c - 1;

% if g(indc) > 0, place this example into the reserve set directly
if (g(indc) >= 0)
   
   % move the example to the reserve set
   bookkeeping(indc,UNLEARNED,RESERVE);
   nstatus = RESERVE;
   
   return;
end;

% compute Qcc and Qc if necessary
num_MVs = length(ind{MARGIN});
Qc = cell(3,1);
if (num_MVs == 0)
	if (length(ind{ERROR}) > 0)
   	Qc{ERROR} = (y(ind{ERROR})*y(indc)).*kernel(X(:,ind{ERROR}),X(:,indc),type,scale);
   end;
else
	Qc{MARGIN} = (y(ind{MARGIN})*y(indc)).*K(1:num_MVs);
	if (length(ind{ERROR}) > 0)
   	Qc{ERROR} = (y(ind{ERROR})*y(indc)).*K(num_MVs+1:length(K));
	end;
end;
if (length(ind{RESERVE}) > 0)
   Qc{RESERVE} = (y(ind{RESERVE})*y(indc)).*kernel(X(:,ind{RESERVE}),X(:,indc),type,scale);
end;
Qcc = kernel(X(:,indc),X(:,indc),type,scale) + deps;

converged = 0;
while (~converged)
   
   perturbations = perturbations + 1;
   
   if (num_MVs > 0)  % change in alpha_c permitted
   
      % compute Qc, beta and gamma
      beta = -Rs*[y(indc) ; Qc{MARGIN}];
      gamma = zeros(size(Q,2),1);
      ind_temp = [ind{ERROR} ind{RESERVE} indc];
      gamma(ind_temp) = [Qc{ERROR} ; Qc{RESERVE} ; Qcc] + Q(:,ind_temp)'*beta;
      
      % check if gamma_c < 0 (kernel matrix is not positive semi-definite)
      if (gamma(indc) < 0)
         error('LEARN: gamma_c < 0');
      end;
      
   else  % change in alpha_c not permitted since the constraint on the sum of the
         % alphas must be preserved.  only b can change.  
      
      % set beta and gamma
      beta = y(indc);
      gamma = y(indc)*y;
      
   end;
   
   % minimum acceptable parameter change (change in alpha_c (num_MVs > 0) or b (num_MVs = 0))
   [min_delta_param,indss,cstatus,nstatus] = min_delta_acb(indc,gamma,beta,1,rflag);
   
   % update a, b, and g
   if (num_MVs > 0)
      a(indc) = a(indc) + min_delta_param;
      a(ind{MARGIN}) = a(ind{MARGIN}) + beta(2:num_MVs+1)*min_delta_param;
   end;   
   b = b + beta(1)*min_delta_param;
   g = g + gamma*min_delta_param;
         
   % update Qc and perform bookkeeping         
   converged = (indss == indc);
   if (converged)
      cstatus = UNLEARNED;
     	Qc{nstatus} = [Qc{nstatus} ; Qcc];
  	else
  		ind_temp = find(ind{cstatus} == indss);
  		Qc{nstatus} = [Qc{nstatus} ; Qc{cstatus}(ind_temp)];
  		Qc{cstatus}(ind_temp) = [];
   end;
   [indco,removed_i] = bookkeeping(indss,cstatus,nstatus);
   if ((nstatus == RESERVE) & (removed_i > 0))
      Qc{nstatus}(removed_i) = [];
   end;
      
   % set g(ind{MARGIN}) to zero
   g(ind{MARGIN}) = 0;
   
   % update Rs and Q if necessary
   if (nstatus == MARGIN)
              
      num_MVs = num_MVs + 1;
      if (num_MVs > 1)
         if (converged)
            gamma = gamma(indss);
         else
               
            % compute beta and gamma for indss            
            beta = -Rs*Q(:,indss);
            gamma = kernel(X(:,indss),X(:,indss),type,scale) + deps + Q(:,indss)'*beta;
            
         end;
      end;
            
      % expand Rs and Q
      updateRQ(beta,gamma,indss);
      
   elseif (cstatus == MARGIN)      
              
      % compress Rs and Q      
      num_MVs = num_MVs - 1;
      updateRQ(indco);
            
   end;         
   
end;

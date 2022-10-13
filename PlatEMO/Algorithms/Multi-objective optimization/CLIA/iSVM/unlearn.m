% UNLEARN - Unlearns the given example if it is a margin or error vector. 
%           If the example is a reserve vector, the SVM remains unchanged
%           and the example status is changed to unlearned.  
%
% Syntax: trans = unlearn(indc)
%
%   indc: index of the example to decrement
%  trans: list of example transitions (index,current status,new status)
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function trans = unlearn(indc)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

% define global variables 
global a;      % alpha coefficients
global b;      % bias
global C;      % regularization parameters
global deps;   % jitter factor in kernel matrix
global g;      % partial derivatives of cost function w.r.t. alpha coefficients
global ind;    % cell array containing indices of margin, error, reserve and unlearned vectors
global Q;      % extended kernel matrix for all vectors
global Rs;     % inverse of extended kernel matrix for margin vectors   
global scale;  % kernel scale
global type;   % kernel type
global X;      % matrix of margin, error, reserve and unlearned vectors stored columnwise
global y;      % column vector of class labels (-1/+1) for margin, error, reserve and unlearned vectors

trans = [];

num_MVs = length(ind{MARGIN});
if (g(indc) < 0)
   
   % remove indc from the indices of error vectors
   i = find(indc == ind{ERROR});
   ind{ERROR}(i) = [];
   
elseif (g(indc) == 0)
   
   % check to see whether indc is labeled as a margin or reserve vector
   % (when learning examples, if g >= 0, they are immediately labeled
   % reserve vectors)
   i = find(indc == ind{MARGIN});
   ismargin = (length(i) > 0);
   
   if (ismargin)
      ind{MARGIN}(i) = [];
      num_MVs = num_MVs - 1;
   
      % do not enforce the margin vector constraint for the example
      % being unlearned
      updateRQ(i+1);
  else
      i = find(indc == ind{RESERVE});
      ind{RESERVE}(i) = [];
  end;
   
else
   
   % remove indc from the indices of reserve vectors
   i = find(indc == ind{RESERVE});
   ind{RESERVE}(i) = [];
   
   % add indc to the indices of unlearned vectors
   ind{UNLEARNED} = [ind{UNLEARNED} indc];
   
end;

if (g(indc) <= 0)   
   
   % compute Qcc and Qc if necessary
   Qc = cell(4,1);
   if (num_MVs > 0)
 	   Qc{MARGIN} = (y(ind{MARGIN})*y(indc)).*kernel(X(:,ind{MARGIN}),X(:,indc),type,scale);
   end;
   if (length(ind{ERROR}) > 0)
  	   Qc{ERROR} = (y(ind{ERROR})*y(indc)).*kernel(X(:,ind{ERROR}),X(:,indc),type,scale);
   end;
   if (length(ind{RESERVE}) > 0)
      Qc{RESERVE} = (y(ind{RESERVE})*y(indc)).*kernel(X(:,ind{RESERVE}),X(:,indc),type,scale);
   end;
   if (length(ind{UNLEARNED}) > 0)
      Qc{UNLEARNED} = (y(ind{UNLEARNED})*y(indc)).*kernel(X(:,ind{UNLEARNED}),X(:,indc),type,scale);
   end;
   Qcc = kernel(X(:,indc),X(:,indc),type,scale) + deps;

   % unlearn indc
   converged = 0;
   while (~converged)
         
      if (num_MVs > 0)  % change in alpha_c permitted
   
         % compute Qc, beta and gamma
		 beta = -Rs*[y(indc) ; Qc{MARGIN}];
         gamma = zeros(size(Q,2),1);
         ind_temp = [ind{ERROR} ind{RESERVE} ind{UNLEARNED} indc];
		 gamma(ind_temp) = [Qc{ERROR} ; Qc{RESERVE} ; Qc{UNLEARNED} ; Qcc] + Q(:,ind_temp)'*beta;
      
      else  % change in alpha_c not permitted since the constraint on the sum of the
            % alphas must be preserved.  only b can change.  
         
         % set beta and gamma
         beta = y(indc);
         gamma = y(indc)*y;
         
      end;
      
      % minimum acceptable parameter change (change in alpha_c (num_MVs > 0) or b (num_MVs = 0))
      [min_delta_param,indss,cstatus,nstatus] = min_delta_acb(indc,gamma,beta,-1,1);
      converged = (indss == indc);
      trans = [trans ; [indss cstatus nstatus]];
      
      % update a, b and g
      if (num_MVs > 0)
         a(indc) = a(indc) + min_delta_param;
         a(ind{MARGIN}) = a(ind{MARGIN}) + beta(2:num_MVs+1)*min_delta_param;
      end;   
      b = b + beta(1)*min_delta_param;
      g = g + gamma*min_delta_param;
   
      if (~converged)
         
	      % update Qc and perform bookkeeping         
  		   ind_temp = find(ind{cstatus} == indss);
		   Qc{nstatus} = [Qc{nstatus} ; Qc{cstatus}(ind_temp)];
 		   Qc{cstatus}(ind_temp) = [];
   	   [indco,removed_i] = bookkeeping(indss,cstatus,nstatus);
         if ((nstatus == RESERVE) & (removed_i > 0))
            Qc{nstatus}(removed_i) = [];
         end;
      
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
      else
      
         % add indc to the appropriate list of indices
         ind{nstatus} = [ind{nstatus} indc];      
      
         if (nstatus == MARGIN)
         
            num_MVs = num_MVs + 1;
            if (num_MVs > 1)
         
               % compute beta and gamma for indss            
               beta = -Rs*Q(:,indss);
               gamma = kernel(X(:,indss),X(:,indss),type,scale) + deps + Q(:,indss)'*beta;
            
            end;
          
            % expand Rs and Q
            updateRQ(beta,gamma,indss);
         
         end;
      end;
   
      % set g(ind{MARGIN}) to zero
      g(ind{MARGIN}) = 0;
   
   end;
end;

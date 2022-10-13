% GESTLOO - Estimates g for the given example that results after 
%           unlearning the example.
%
% Syntax: g_est = gestloo(indc)
%
%   indc: example index
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function g_est = gestloo(indc)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

% define global variables 
global a;      % alpha coefficients
global C; 	   % regularization parameters
global deps;   % jitter factor in kernel matrix
global g;      % partial derivatives of cost function w.r.t. alpha coefficients
global ind;    % cell array containing indices of margin, error, reserve and unlearned vectors
global Q;      % extended kernel matrix for all vectors
global Rs;     % inverse of extended kernel matrix for margin vectors   
global scale;  % kernel scale
global type;   % kernel type
global X;      % matrix of margin, error, reserve and unlearned vectors stored columnwise

num_MVs = length(ind{MARGIN});
i = find(indc == ind{MARGIN});
is_margin = ~isempty(i);
if (is_margin)
       
   % shrink the inverse of the extended kernel matrix for the margin vectors.
   % assumes there is more than one margin vector.
   stripped = [1:i i+2:size(Rs,1)]; 
   Rs_new = Rs(stripped,stripped)-Rs(stripped,i+1)*Rs(i+1,stripped)/Rs(i+1,i+1);
   
   Qcc = Q(i+1,indc);
   Qc = [Q(1:i,indc) ; Q(i+2:size(Q,1),indc)];
        
else
      
   Rs_new = Rs;
   Qcc = kernel(X(:,indc),X(:,indc),type,scale) + deps;
   Qc = Q(:,indc);
   
end;
   
% compute beta and gamma_c
beta = -Rs_new*Qc;
gamma_c = Qcc + Qc'*beta;

% compute g_est
g_est = g(indc) - a(indc)*gamma_c;

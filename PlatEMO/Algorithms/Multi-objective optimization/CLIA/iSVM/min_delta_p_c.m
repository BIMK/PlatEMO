% MIN_DELTA_P_C - Computes the minimum acceptable change in p_c during
%                 regularization parameter perturbation and indicates 
%                 which example changes status along with the type of status 
%                 change.
%
% Syntax: [min_dpc,indss,cstatus,nstatus] = min_delta_p_c(p_c,gamma,beta,lambda)
%
%    min_dpc: minimum acceptable change in p_c
%      indss: the example changing status (0 if nothing changed status)
%    cstatus: current status of the example
%             0: if indss = indc
%             1: margin vector
%             2: error vector
%             3: reserve vector
%             4: unlearned vector
%    nstatus: new status of the example
%             returns one of the values above
%        p_c: current value of p_c
%      gamma: margin sensitivities (delta g / delta p_c)
%       beta: coefficient sensitivities (delta alpha_k/b / delta p_c)
%     lambda: regularization sensitivities (delta C_k / delta p_c)
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [min_dpc,indss,cstatus,nstatus] = min_delta_p_c(p_c,gamma,beta,lambda)

% flags for example state  
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

% define global variables
global a;       % the alpha coefficients
global C;       % regularization parameter(s)
global g;       % partial derivatives of the objective function w.r.t. the alphas
global ind;     % cell array containing indices of margin, error, reserve and unlearned vectors

indss = zeros(5,1);
cstatus = zeros(5,1);
nstatus = zeros(5,1);

% upper limit on change in p_c assuming no other examples change status 
delta_p_c = 1-p_c;

% change in p_c that causes a margin vector to change to a reserve vector
if (length(beta) > 1)  % if there are margin vectors  
   beta_s = beta(2:length(beta));
   flags = (beta_s < 0);
   [delta_mr,i] = min_delta(flags,a(ind{MARGIN}),zeros(size(a(ind{MARGIN}))),beta_s);
   if (delta_mr < Inf)
      indss(2) = ind{MARGIN}(i); 
      cstatus(2) = MARGIN;
      nstatus(2) = RESERVE;   
   end;
else
   delta_mr = Inf;
end;

% change in p_c that causes a margin vector to change to an error vector
if (length(beta) > 1)  % if there are margin vectors
   lambda_s = lambda(ind{MARGIN});
   v = beta_s-lambda_s;
   flags = (v > eps);
   if (any(flags))
      not_z = find(v > 0);
      delta_me = Inf*ones(size(v));
      delta_me(not_z) = (C(ind{MARGIN}(not_z))-a(ind{MARGIN}(not_z)))./v(not_z);
      [delta_me,i] = min(delta_me);
      if (delta_me < Inf)
         indss(3) = ind{MARGIN}(i);
         cstatus(3) = MARGIN;
         nstatus(3) = ERROR;
      end;
   else
      delta_me = Inf;
   end;
else
   delta_me = Inf;   
end;

% change in p_c that causes an error vector to change to a margin vector
gamma_e = gamma(ind{ERROR});
flags = (gamma_e > 0);
[delta_em,i] = min_delta(flags,g(ind{ERROR}),zeros(length(ind{ERROR}),1),gamma_e);
if (delta_em < Inf)
   indss(4) = ind{ERROR}(i);
   cstatus(4) = ERROR;
   nstatus(4) = MARGIN;
end;

% change in p_c that causes a reserve vector to change to a margin vector
gamma_r = gamma(ind{RESERVE});
flags = (g(ind{RESERVE}) >= 0) & (gamma_r < 0);
[delta_rm,i] = min_delta(flags,g(ind{RESERVE}),zeros(length(ind{RESERVE}),1),gamma_r);
if (delta_rm < Inf)
   indss(5) = ind{RESERVE}(i);
   cstatus(5) = RESERVE;
   nstatus(5) = MARGIN;    
end;

% minimum acceptable value for p_c
[min_dpc,min_ind] = min([delta_p_c,delta_mr,delta_me,delta_em,delta_rm]);
indss = indss(min_ind);
cstatus = cstatus(min_ind);
nstatus = nstatus(min_ind);


% MIN_DELTA_P_S - Computes the minimum acceptable change in p_s during
%                 kernel parameter perturbation and indicates which 
%                 example changes status along with the type of status change.
%
% Syntax: [min_dps,indss,cstatus,nstatus] = min_delta_p_s(p_s,gamma,beta,lambda)
%
%    min_dps: minimum acceptable change in p_s
%      indss: the example changing status (0 if nothing changed status)
%    cstatus: current status of the example
%             0: if indss = indc
%             1: margin vector
%             2: error vector
%             3: reserve vector
%             4: unlearned vector
%    nstatus: new status of the example
%             returns one of the values above
%        p_s: current value of p_s
%      gamma: margin sensitivities (delta g / delta p_s)
%       beta: coefficient sensitivities (delta alpha_k/b / delta p_s)
%     lambda: regularization sensitivities (delta C_k / delta p_s)
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [min_dps,indss,cstatus,nstatus] = min_delta_p_s(p_s,gamma,beta,lambda)

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

indss = zeros(6,1);
cstatus = zeros(6,1);
nstatus = zeros(6,1);

if (length(ind{UNLEARNED}) > 0)
   
   delta_p_s = Inf;
   
   % change in p_s that causes an unlearned vector to change to an error vector
   lambda_u = lambda(ind{UNLEARNED});
   flags = (lambda_u > 0);
   [delta_ue,i] = min_delta(flags,a(ind{UNLEARNED}),C(ind{UNLEARNED}),lambda_u);
   if (delta_ue < Inf)
      indss(1) = ind{UNLEARNED}(i); 
      cstatus(1) = UNLEARNED;
      nstatus(1) = ERROR;   
   end;
   
   % change in p_s that causes an unlearned vector to change to a margin vector
   gamma_u = gamma(ind{UNLEARNED});
   g_u = g(ind{UNLEARNED});
   flags = (sign(gamma_u).*sign(g_u) == -1);
   [delta_um,i] = min_delta(flags,g_u,zeros(size(g_u)),gamma_u);
   if (delta_um < Inf)
      indss(2) = ind{UNLEARNED}(i);
      cstatus(2) = UNLEARNED;
      nstatus(2) = MARGIN;
   end;
      
   % change in p_s that causes an unlearned vector to change to a reserve vector
   flags = (lambda_u < 0);
   [delta_ur,i] = min_delta(flags,a(ind{UNLEARNED}),zeros(length(ind{UNLEARNED}),1),lambda_u);
   if (delta_ur < Inf)
      indss(3) = ind{UNLEARNED}(i); 
      cstatus(3) = UNLEARNED;
      nstatus(3) = RESERVE;   
   end;
   
else
   delta_p_s  = 1-p_s;
   delta_ue = Inf;
   delta_um = Inf;
   delta_ur = Inf;
end;

% change in p_s that causes a margin vector to change to an error vector
% or reserve vector
if (length(beta) > 1)  % if there are margin vectors  
   beta_s = beta(2:length(beta));
   flags = (abs(beta_s) > 0);
   [delta_mer,i] = min_delta(flags,a(ind{MARGIN}),C(ind{MARGIN}).*(beta_s > 0),beta_s);
   if (delta_mer < Inf)
      indss(5) = ind{MARGIN}(i); 
      cstatus(5) = MARGIN;
      nstatus(5) = ERROR*(beta_s(i) > 0) + RESERVE*(beta_s(i) < 0);   
   end;
else
   delta_mer = Inf;
end;

% change in p_s that causes an error vector to change to a margin vector
gamma_e = gamma(ind{ERROR});
flags = (gamma_e > 0);
[delta_em,i] = min_delta(flags,g(ind{ERROR}),zeros(length(ind{ERROR}),1),gamma_e);
if (delta_em < Inf)
   indss(6) = ind{ERROR}(i);
   cstatus(6) = ERROR;
   nstatus(6) = MARGIN;
end;

% change in p_s that causes a reserve vector to change to a margin vector
gamma_r = gamma(ind{RESERVE});
flags = (g(ind{RESERVE}) >= 0) & (gamma_r < 0);
[delta_rm,i] = min_delta(flags,g(ind{RESERVE}),zeros(length(ind{RESERVE}),1),gamma_r);
if (delta_rm < Inf)
   indss(7) = ind{RESERVE}(i);
   cstatus(7) = RESERVE;
   nstatus(7) = MARGIN;    
end;

% minimum acceptable value for p_s
[min_dps,min_ind] = min([delta_ue,delta_um,delta_ur,delta_p_s,delta_mer,delta_em,delta_rm]);
indss = indss(min_ind);
cstatus = cstatus(min_ind);
nstatus = nstatus(min_ind);


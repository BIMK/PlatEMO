% MIN_DELTA_ACB - Computes the minimum acceptable change in alpha_c or b
%                 and indicates which example changes status along
%                 with the type of status change.
%
% Syntax: [min_dacb,indss,cstatus,nstatus] = min_delta_acb(indc,gamma,beta,polc,rflag)
%
%   min_dacb: minimum acceptable change in alpha_c (num MVs > 0) / b (num MVs = 0)
%      indss: the example changing status
%    cstatus: current status of the example
%             0: if indss = indc
%             1: margin vector
%             2: error vector
%             3: reserve vector
%             4: unlearned vector
%    nstatus: new status of the example
%             returns one of the values above
%       indc: the current example being learned/unlearned
%      gamma: margin sensitivities (delta g / delta alpha_c/b)
%       beta: coefficient sensitivities (delta alpha_k/b / delta alpha_c/b)
%       polc: sign of the change in alpha_c/b (+1/-1)
%      rflag: flag indicating whether or not to check if any reserve vectors
%             become margin vectors during learning
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [min_dacb,indss,cstatus,nstatus] = min_delta_acb(indc,gamma,beta,polc,rflag)

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
indss(1:3) = indc;
cstatus = zeros(6,1);
nstatus = zeros(6,1);

% upper limits on change in alpha_c or b assuming no other examples change status 
if (polc == 1)  
   delta_m = -g(indc)/gamma(indc);
   nstatus(1) = MARGIN;   
   if (length(beta) > 1)   % if there are margin vectors 
      delta_e = C(indc) - a(indc);
      nstatus(2) = ERROR;   
   else                    % only the bias term is allowed to change so indc can't be an error vector  
      delta_e = Inf;
   end;
   delta_r = Inf;
else  
   if (g(indc) > 0)        % decrementing indc when a(indc) > 0 and g(indc) > 0 during kernel perturbation
      delta_m = g(indc)/gamma(indc);
      nstatus(1) = MARGIN;
   else
      delta_m = Inf;
   end;
   if (a(indc) <= C(indc))
      delta_e = Inf;
   	  delta_r = a(indc);
      if (g(indc) > 0)     % decrementing indc when a(indc) > 0 and g(indc) > 0 during kernel perturbation
         nstatus(3) = RESERVE;
      else
         nstatus(3) = UNLEARNED;
      end;
   else                    % decrementing indc when a(indc) > C(indc) during reg. parameter perturbation
      delta_e = a(indc) - C(indc);
      delta_r = Inf;
      nstatus(2) = ERROR;
   end;
end;

% change in alpha_c or b that causes a margin vector to change to an error vector
% or reserve vector
if (length(beta) > 1)  % if there are margin vectors  
   beta_s = polc*beta(2:length(beta));
   flags = (abs(beta_s) > 0);
   [delta_mer,i] = min_delta(flags,a(ind{MARGIN}),C(ind{MARGIN}).*(beta_s > 0),beta_s);
   if (delta_mer < Inf)
      indss(4) = ind{MARGIN}(i); 
      cstatus(4) = MARGIN;
      nstatus(4) = ERROR*(beta_s(i) > 0) + RESERVE*(beta_s(i) < 0);   
   end;
else
   delta_mer = Inf;
end;

% change in alpha_c or b that causes an error vector to change to a margin vector
gamma_e = polc*gamma(ind{ERROR});
flags = (gamma_e > 0);
[delta_em,i] = min_delta(flags,g(ind{ERROR}),zeros(length(ind{ERROR}),1),gamma_e);
if (delta_em < Inf)
   indss(5) = ind{ERROR}(i);
   cstatus(5) = ERROR;
   nstatus(5) = MARGIN;
end;

% change in alpha_c or b that causes a reserve vector to change to a margin vector
if (rflag)
   gamma_r = polc*gamma(ind{RESERVE});
   flags = (g(ind{RESERVE}) >= 0) & (gamma_r < 0);
   [delta_rm,i] = min_delta(flags,g(ind{RESERVE}),zeros(length(ind{RESERVE}),1),gamma_r);
   if (delta_rm < Inf)
      indss(6) = ind{RESERVE}(i);
      cstatus(6) = RESERVE;
      nstatus(6) = MARGIN;    
   end;
else
   delta_rm = Inf;
end;

% minimum acceptable value for |delta_ac| or |delta_b|
[min_dacb,min_ind] = min([delta_m,delta_e,delta_r,delta_mer,delta_em,delta_rm]);
indss = indss(min_ind);
cstatus = cstatus(min_ind);
nstatus = nstatus(min_ind);

% multiply by the proper sign to yield delta_ac or delta_b
min_dacb = polc*min_dacb;

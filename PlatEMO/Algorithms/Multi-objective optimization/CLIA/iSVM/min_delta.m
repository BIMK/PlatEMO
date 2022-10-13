% MIN_DELTA - Computes the minimum change in the parameter lambda
%             that causes one of the parameters psi(k) to change from
%             the given initial value to the given final value.
%             Only the parameters with flags(k) = 1 are checked.  
%             psi(k) and lambda are assumed to be linearly related.  
%
% Syntax: [min_d,k] = min_delta(flags,psi_initial,psi_final,psi_sens)
%
%       flags: parameters to check
%       min_d: minimum change in lambda
%           k: parameter psi_k that achieves psi_final_k
% psi_initial: initial values for the parameters psi_i
%   psi_final: final values for the parameters psi_i
%    psi_sens: parameter sensitivities (d(psi_i)/d(lambda))
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [min_d,k] = min_delta(flags,psi_initial,psi_final,psi_sens)

if (any(flags))
   
   % find the parameters to check
   ind = find(flags);

   deltas = (psi_final(ind)-psi_initial(ind))./psi_sens(ind);
   [min_d,i] = min(deltas);
   k = ind(i);

   % if there is more than one parameter that achieves the final value
   % with min_delta, select the parameter with the highest 
   % sensitivity |psi_sens|
   min_k = find(deltas == min_d);
   if (length(min_k) > 1)
      [max_sens,i] = max(abs(psi_sens(ind(min_k))));
      k = ind(min_k(i));
   end;

else
   
   min_d = Inf;
   k = -1;
   
end;

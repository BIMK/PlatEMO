% MOVE_INDR - Removes the indices indc from the indices inda
%             and adds them to the reserve vector index list
%             if necessary.  
%
% Syntax: [inda,i] = move_indr(inda,indc)
%
%      i: when i > 0, the i-th element was removed from ind{RESERVE}
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [inda,removed_i] = move_indr(inda,indc)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLABELED = 4;

% define global variables
global g;							% partial derivatives of cost function w.r.t. alpha coefficients
global ind;							% cell array containing indices of margin, error, reserve and unlearned vectors
global max_reserve_vectors;   % maximum number of reserve vectors

removed_i = 0;
num_RVs_orig = length(ind{RESERVE});

% shift indc from inda to ind{RESERVE}
[inda,ind{RESERVE}] = move_ind(inda,ind{RESERVE},indc);

% if we need to remove some reserve vectors
if (length(ind{RESERVE}) > max_reserve_vectors) 
   
   % sort g(ind{RESERVE})
   [g_sorted,i] = sort(g(ind{RESERVE}));
   
   % reserve vectors that need to be removed
   removed = i(max_reserve_vectors+1:length(i));
   
   % find any original reserve vectors that need to be removed
   k = find(removed <= num_RVs_orig);
   if (length(k) > 0)
      removed_i = removed(k);
   end;
   
   % remove the necessary reserve vectors
   ind{RESERVE}(removed) = [];
   
end;


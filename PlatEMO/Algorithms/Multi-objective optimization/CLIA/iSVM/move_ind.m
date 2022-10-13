% MOVE_IND - Removes the indices indc from the indices inda
%            and appends them to end of the indices indb.
%            The relative order of the remaining indices in 
%            inda is preserved.  
%
% Syntax: [inda,indb] = move_ind(inda,indb,indc)
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function [inda,indb] = move_ind(inda,indb,indc)

if (~isempty(indc))
   indb = [indb indc];
   new_inds = [];
   for i = 1:length(inda)
      if (~any(inda(i) == indc))
         new_inds = [new_inds inda(i)];
      end;
   end;
   inda = new_inds;
end;

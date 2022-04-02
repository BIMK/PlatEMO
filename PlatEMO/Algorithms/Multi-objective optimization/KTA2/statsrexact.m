function [pval,P] = statsrexact(v,w)
%STATSREXACT Compute exact tail probability for signed rank statistic.
%   [PVAL,ALLP]=STATSREXACT(V,W) computes the tail probability PVAL
%   for the statistic W with the vector V of ranks.  ALLP is a matrix
%   containing the probabilities (col. 2) for each W value (col. 1).
%
%   Private function used by the SIGNRANK function.

%   Copyright 2003-2006 The MathWorks, Inc. 


n = length(v);
v = sort(v(:)'); % make sure it's a row

% For convenience we can just compute the lower tail.  If w is
% in the upper tail, compute its equivalent lower tail value.
maxw = n*(n+1)/2;
folded = (w>maxw/2);
if folded
   w = maxw-w;
end

% We would like to use the elements of w and v as indexes into
% arrays that enumerate possible values.  If there are ties causing
% non-integer ranks, multiply by 2 to force everything to integer.
doubled = any(v~=floor(v));
if doubled
   v = round(2*v);
   w = round(2*w);
end

C = zeros(w+1,1);  % C(w+1) will be the number of combinations adding
                   % to w at each step
C(1) = 1;          % just one combination includes nothing
top = 1;           % top entry currently in use in C vector

% Look at all combinations of ranks that could contribute
% to the observed value of W
for vj=v(v<=w)

   % C now enumerates combinations not including v(j).  Now update the
   % elements that could include v(j).
   newtop = min(top+vj,w+1);
   hi = min(vj,w+1)+1:newtop;
   lo = 1:length(hi);

   C(hi) = C(hi) + C(lo);

   top = newtop;
end

% Convert to probabilities
C = C / (2^n);

% Get tail probability
pval = sum(C);

if nargout>1
   allw = 0:w;
   if doubled
      allw = allw/2;
   end
   if folded
      allw = n*(n+1)/2 - allw;
   end
   
   P = [allw(:), C(:)];
end
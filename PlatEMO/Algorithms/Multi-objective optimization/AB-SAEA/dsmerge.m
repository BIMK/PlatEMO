function  [mS, mY] = dsmerge(S, Y, ds, nms, wtds, wtdy)
%DSMERGE  Merge data for multiple design sites.
%
% Call
%   [mS, mY] = dsmerge(S, Y)
%   [mS, mY] = dsmerge(S, Y, ds)
%   [mS, mY] = dsmerge(S, Y, ds, nms)
%   [mS, mY] = dsmerge(S, Y, ds, nms, wtds)
%   [mS, mY] = dsmerge(S, Y, ds, nms, wtds, wtdy)
%
% Input
% S, Y : Data points (S(i,:), Y(i,:)), i = 1,...,m
% ds   : Threshold for equal, normalized sites. Default is 1e-14.
% nms  : Norm, in which the distance is measured.
%        nms =  1 : 1-norm (sum of absolute coordinate differences)
%               2 : 2-norm (Euclidean distance) (default)
%        otherwise: infinity norm (max coordinate difference)      
% wtds : What to do with the S-values in case of multiple points.
%        wtds = 1 : return the mean value (default)
%               2 : return the median value
%               3 : return the 'cluster center'
% wtdy : What to do with the Y-values in case of multiple points.
%        wtdy = 1 : return the mean value (default)
%               2 : return the median value
%               3 : return the 'cluster center' value
%               4 : return the minimum value
%               5 : return the maximum value    
%
% Output
% mS : Compressed design sites, with multiple points merged
%      according to wtds
% mY : Responses, compressed according to wtdy

% hbn@imm.dtu.dk  
% Last update July 3, 2002

% Check design points
[m n] = size(S);  % number of design sites and their dimension
sY = size(Y);
if  min(sY) == 1,  Y = Y(:);   lY = max(sY);  sY = size(Y);
else,              lY = sY(1); end
if m ~= lY
  error('S and Y must have the same number of rows'), end

% Threshold
if  nargin < 3
  ds = 1e-14;
elseif  (ds < 0) | (ds > .5)
  error('ds must be in the range [0, 0.5]'), end

% Which measure
if  nargin < 4
  nms = 2;
elseif  (nms ~= 1) & (nms ~= 2)
  nms = Inf;
end

% What to do
if  nargin < 5
  wtds = 1;
else
  wtds = round(wtds);
  if  (wtds < 1) | (wtds > 3)
    error('wtds must be in the range [1, 3]'), end
end
if  nargin < 6
  wtdy = 1;
else
  wtdy = round(wtdy);
  if  (wtdy < 1) | (wtdy > 5)
    error('wtdy must be in the range [1, 5]'), end
end

% Process data
more = 1;
ladr = zeros(1,ceil(m/2));
while more
  m = size(S,1);
  D = zeros(m,m);
  
  % Normalize sites
  mS = mean(S);   sS = std(S);
  scS = (S - repmat(mS,m,1)) ./ repmat(sS,m,1);
  
  % Calculate distances D (upper triangle of the symetric matrix)
  for k = 1 : m-1
    kk = k+1 : m;
    dk = abs(repmat(scS(k,:), m-k, 1) - scS(kk,:));
    if  nms == 1,      D(kk,k) = sum(dk,2);
    elseif  nms == 2,  D(kk,k) = sqrt(sum(dk.^2,2));
    else,              D(kk,k) = max(dk,[],2);      end
  end
  D = D + D'; % make D symetric
  
  % Check distances
  mult = zeros(1,m);
  for  j = 1 : m
    % Find the number of multiple sites in each column of D
    mult(j) = length(find(D(:,j) < ds));
  end
  % Find the first column with the maximum number of multiple sites
  [mmult jj] = max(mult);
  
  if  mmult == 1,  more = 0;
  else
    nm = 0;
    while  mmult > 1
      nm = nm + 1;  % no. of points to merge
      ladr(nm) = jj;
      
      % Merge point no jj and its neighbours, note that jj is the center
      % of the cluster, as it has the most neighbors (among the multiple sites)
      ngb = find(D(:,jj) < ds);
      
      switch  wtds
      case 1,  S(jj,:) = mean(S(ngb,:));
      case 2,  S(jj,:) = median(S(ngb,:));
      case 3,  S(jj,:) = S(jj,:);
      end
      
      switch  wtdy
      case 1,  Y(jj,:) = mean(Y(ngb,:));
      case 2,  Y(jj,:) = median(Y(ngb,:));
      case 3,  Y(jj,:) = Y(jj,:);
      case 4,  Y(jj,:) = min(Y(ngb,:));
      case 5,  Y(jj,:) = max(Y(ngb,:));
      end
      
      % Delete from list
      mult(ngb) = 0;
      [mmult jj] = max(mult);
    end

    % Reduced data set
    act = [find(mult > 0)  ladr(1:nm)];
    S = S(act,:);    Y = Y(act,:);
  end % multiple
end % loop

% Return reduced set
mS = S;   mY = Y;
function [I,I2] = clusterbal(data, ndata, k, ncenters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

niter = 300;

[centers, idxs] = anderbergy(data, ncenters);

id   = eye(ncenters);
ea   = 1000000; 
maxi = 1000000;
for kk = 1 : niter
    d2  = dist2(data, centers);  
    d2q = d2;
    
    % Assign each point to nearest centre
    [minvals, index] = min(d2', [], 1);
    post = id(index,:);
    nat  = 1 : ndata;
    aux  = [nat',minvals'];
    aux1 = sortrows(aux,2);
    idx  = aux1(:,1);	
  
    asign = k*ones(1,ncenters);
    postf = zeros(ndata,ncenters);
    for j = 1 : ndata;
        
      i = idx(j);
      [minrow, icol] = min (d2q(i,:));
      asign(icol)    = asign(icol) - 1;
      postf(i,icol)  = 1;
      
      if asign(icol) == 0 
          d2q(:,icol) = maxi;
      end
      
    end

    post = postf;
    minima = d2.*post;
    e = sum( sum (minima) ); 
    if e >= ea
        post    = posta;
        e       = ea;
        d2      = d2a;
        centers = centersa;
    else
        ea;
        kk;
    end
    if kk==niter 
        a = [data,post];
        break
    end

    ea       = e;
    posta    = post;
    d2a      = d2;
    centersa = centers;

    [maxc,I,J] = maxmatrix(d2);
    auxi = randperm(ncenters);
    if auxi(1)~=J
      H = auxi(1);
    else
      H = auxi(2);
    end
    
    [maxaux,L] = max(d2(:,H));
    post(I,J)  = 0;
    post(I,H)  = 1;
    post(L,H)  = 0;
    post(L,J)  = 1;

    for j = 1 : ncenters
        centers(j,:) = sum( data( find( post(:,j) ), : ), 1)/k;
    end  
end

I = [];
for i = 1:ncenters
    I=[I find(post(:,i))'];
end

I2 = [];
for i = 1:ndata
    I2 = [I2 find(post(i,:))];
end

I  = I.';
I2 = I2.';

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n2 = dist2(x, c)
%DIST2	Calculates squared distance between two sets of points.
%
%	Description
%	D = DIST2(X, C) takes two matrices of vectors and calculates the
%	squared Euclidean distance between them.  Both matrices must be of
%	the same column dimension.  If X has M rows and N columns, and C has
%	L rows and N columns, then the result has M rows and L columns.  The
%	I, Jth entry is the  squared distance from the Ith row of X to the
%	Jth row of C.
%
%	See also
%	GMMACTIV, KMEANS, RBFFWD
%

%	Copyright (c) Ian T Nabney (1996-2001)

[ndata, dimx] = size(x);
[ncentres, dimc] = size(c);
if dimx ~= dimc
	error('Data dimension does not match dimension of centres')
end

n2 = (ones(ncentres, 1) * sum((x.^2)', 1))' + ...
  ones(ndata, 1) * sum((c.^2)',1) - ...
  2.*(x*(c'));

% Rounding errors occasionally cause negative entries in n2
if any(any(n2<0))
  n2(n2<0) = 0;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxc, I, J] = maxmatrix(A)

[maxr,I] = max(A);
[maxc,J] = max(maxr);
I = I(J);

return

function [model, time] = rbf_build(Xtr, Ytr, bf_type, bf_c, usePolyPart, verbose)
% RBFBUILD
% Builds a Radial Basis Function (RBF) interpolant using training data
%
% Call
%   [model, time] = rbf_build(Xtr, Ytr, bf_type, bf_c, usePolyPart, verbose)
%   [model, time] = rbf_build(Xtr, Ytr, bf_type, bf_c, usePolyPart)
%   [model, time] = rbf_build(Xtr, Ytr, bf_type, bf_c)
%   [model, time] = rbf_build(Xtr, Ytr, bf_type)
%   [model, time] = rbf_build(Xtr, Ytr)
%
% Input
% Xtr, Ytr    : Training data points (Xtr(i,:), Ytr(i)), i = 1,...,n
%               Note that the input variables must be scaled to e.g. [0,1]
%               or [-1,1] for better predictive performance.
% bf_type     : Type of the basis functions (default = 'MQ'):
%               'BH' = Biharmonic
%               'MQ' = Multiquadric
%               'IMQ' = Inverse Multiquadric
%               'TPS' = Thin plate spline
%               'G' = Gaussian
% bf_c        : Parameter c value (default = 1)
% usePolyPart : Use also the polynomial term P of the model Y = P + RBF
%               (default = 0, do not use)
% verbose     : Set to 0 for no verbose (default = 1)
%
% Output
% model     : RBF model - a struct with the following elements:
%    n      : Number of data points in the training data set
%    meanY  : Mean of Ytr
%    bf_type: Type of the basis functions
%    bf_c   : Parameter c value
%    poly   : Use also the polynomial term
%    coefs  : Coefficients of the model
% time      : Execution time
%
% Please give a reference to the software web page in any publication
% describing research performed using the software, e.g. like this:
% Jekabsons G. Radial Basis Function interpolation for Matlab, 2009,
% available at http://www.cs.rtu.lv/jekabsons/

% This source code is tested with Matlab version 7.1 (R14SP3).

% =========================================================================
% RBF interpolation
% Version: 1.1
% Date: August 12, 2009
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009  Gints Jekabsons
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% =========================================================================

if nargin < 2
    error('Too few input arguments.');
else

    [n, d] = size(Xtr);
    [ny, dy] = size(Ytr);
    if (n < 2) || (d < 1) || (ny ~= n) || (dy ~= 1)
        error('Wrong training data sizes.');
    end

    if nargin < 3
        bf_type = 'MQ';
    end
    if nargin < 4
        bf_c = 1;
    end
    if nargin < 5
        usePolyPart = 1;
    end
    if nargin < 6
        verbose = 0;
    end

    tic;

    model.n = n;
    model.meanY = mean(Ytr);
    model.bf_type = bf_type;
    model.bf_c = bf_c;
    model.poly = usePolyPart;

    %calculate and transform distances between all the points in the training data
    dist = zeros(n, n);
    switch upper(model.bf_type)
        case 'BH'
            if verbose
                fprintf('Building RBF (biharmonic) model...\n');
            end
            for i = 1 : n
                %for j = i : n
                %    dist(i, j) = norm(Xtr(i,:) - Xtr(j,:));
                %end
                dist(i, i:n) = sqrt(sum((repmat(Xtr(i,:),n-i+1,1) - Xtr(i:n,:)).^2,2));
                dist(i+1:n, i) = dist(i, i+1:n);
            end
        case 'IMQ'
            if verbose
                fprintf('Building RBF (inverse multiquadric) model...\n');
            end
            for i = 1 : n
                %for j = i : n
                %    dist(i, j) = 1 / sqrt(sum((Xtr(i,:) - Xtr(j,:)).^2) + bf_c^2);
                %end
                dist(i, i:n) = 1 ./ sqrt(sum((repmat(Xtr(i,:),n-i+1,1) - Xtr(i:n,:)).^2,2) + bf_c^2);
                dist(i+1:n, i) = dist(i, i+1:n);
            end
        case 'TPS'
            if verbose
                fprintf('Building RBF (thin plate spline) model...\n');
            end
            for i = 1 : n
                %for j = i : n
                %    dist(i, j) = sum((Xtr(i,:) - Xtr(j,:)).^2);
                %    dist(i, j) = (dist(i, j) + bf_c^2) * log(sqrt(dist(i, j) + bf_c^2));
                %end
                dist(i, i:n) = sum((repmat(Xtr(i,:),n-i+1,1) - Xtr(i:n,:)).^2,2);
                dist(i, i:n) = (dist(i, i:n) + bf_c^2) .* log(sqrt(dist(i, i:n) + bf_c^2));
                dist(i+1:n, i) = dist(i, i+1:n);
            end
        case 'G'
            if verbose
                fprintf('Building RBF (Gaussian) model...\n');
            end
            for i = 1 : n
                %for j = i : n
                %    dist(i, j) = exp(-sum((Xtr(i,:) - Xtr(j,:)).^2) / (2*bf_c^2));
                %end
                dist(i, i:n) = exp(-sum((repmat(Xtr(i,:),n-i+1,1) - Xtr(i:n,:)).^2,2) / (2*bf_c^2));
                dist(i+1:n, i) = dist(i, i+1:n);
            end
        otherwise %MQ
            if verbose
                fprintf('Building RBF (multiquadric) model...\n');
            end
            for i = 1 : n
                %for j = i : n
                %    dist(i, j) = sqrt(sum((Xtr(i,:) - Xtr(j,:)).^2) + bf_c^2);
                %end
                dist(i, i:n) = sqrt(sum((repmat(Xtr(i,:),n-i+1,1) - Xtr(i:n,:)).^2,2) + bf_c^2);
                dist(i+1:n, i) = dist(i, i+1:n);
            end
    end

    %calculate coefs
    if model.poly == 0
        model.coefs = dist \ (Ytr - model.meanY);
    else
        A = [dist, ones(n,1), Xtr; [ones(n,1), Xtr]', zeros(d+1,d+1)];
        model.coefs  = A \ [Ytr; zeros(d+1,1)];
    end

    time = toc;

    if verbose
        fprintf('Execution time: %0.2f seconds\n', time);
    end

end
return

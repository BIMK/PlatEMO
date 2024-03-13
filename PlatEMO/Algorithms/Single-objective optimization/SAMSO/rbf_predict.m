function [Yq] = rbf_predict(model, Xtr, Xq)
% RBFPREDICT
% Predicts output values for the given query points Xq using a RBF model
%
% Call
%   [Yq] = rbf_predict(model, Xtr, Xq)
%
% Input
% model     : RBF model
% Xtr       : Inputs of the training data (Xtr(i,:)), i = 1,...,n (the same
%             matrix with which the model was built)
% Xq        : Inputs of query data points (Xq(i,:)), i = 1,...,nq
%
% Output
% Yq        : Predicted outputs of query data points (Yq(i)), i = 1,...,nq
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

if nargin < 3
    error('Too few input arguments.');
end

if model.n ~= size(Xtr, 1)
    error('The matrix Xtr should be the same matrix with which the model was built.');
end

nq = size(Xq, 1);
Yq = zeros(nq, 1);

%dist = zeros(model.n,1);
for t = 1 : nq
    switch upper(model.bf_type)
        case 'BH'
            %for j = 1 : model.n
            %    dist(j) = norm(Xq(t,:) - Xtr(j,:));
            %end
            dist = sqrt(sum((repmat(Xq(t,:),model.n,1) - Xtr).^2,2));
        case 'IMQ'
            %for j = 1 : model.n
            %    dist(j) = 1 / sqrt(sum((Xq(t,:) - Xtr(j,:)).^2) + model.bf_c^2);
            %end
            dist = 1 ./ sqrt(sum((repmat(Xq(t,:),model.n,1) - Xtr).^2,2) + model.bf_c^2);
        case 'TPS'
            %for j = 1 : model.n
            %    dist(j) = sum((Xq(t,:) - Xtr(j,:)).^2);
            %    dist(j) = (dist(j) + model.bf_c^2) * log(sqrt(dist(j) + model.bf_c^2));
            %end
            dist = sum((repmat(Xq(t,:),model.n,1) - Xtr).^2,2);
            dist = (dist + model.bf_c^2) .* log(sqrt(dist + model.bf_c^2));
        case 'G'
            %for j = 1 : model.n
            %    dist(j) = exp(-sum((Xq(t,:) - Xtr(j,:)).^2) / (2*model.bf_c^2));
            %end
            dist = exp(-sum((repmat(Xq(t,:),model.n,1) - Xtr).^2,2) / (2*model.bf_c^2));
        otherwise %MQ
            %for j = 1 : model.n
            %    dist(j) = sqrt(sum((Xq(t,:) - Xtr(j,:)).^2) + model.bf_c^2);
            %end
            dist = sqrt(sum((repmat(Xq(t,:),model.n,1) - Xtr).^2,2) + model.bf_c^2);
    end
    if model.poly == 0
        Yq(t) = model.meanY + model.coefs' * dist;
    else
        Yq(t) = model.coefs' * [dist; 1; Xq(t,:)'];
    end
end

return

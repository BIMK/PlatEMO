function  [y,mse] = predictor(x, dmodel)
% -----------------------------------------------------------------------------------
% this file has been modified to predict multiple test points at the same
% time, x could be a matrix with N * n
% where N is the number of test points to be predicted and n is the
% dimension of design variables.
% MATLAB2016b, zhandawei 2017.04.26
% -----------------------------------------------------------------------------------
%PREDICTOR  Predictor for y(x) using the given DACE model.
%
% Call:   [y,mse] = predictor(x, dmodel)
%
% Input
% x      : trial design sites with n dimensions.
%          For mx trial sites x:
%          If mx = 1, then both a row and a column vector is accepted,
%          otherwise, x must be an mx*n matrix with the sites stored
%          rowwise.
% dmodel : Struct with DACE model; see DACEFIT
%
% Output
% y    : predicted response at x.
% mse  : Estimated mean squared error of the predictor;
% hbn@imm.dtu.dk
% Last update August 26, 2002

[m,n] = size(dmodel.S);  % number of design sites and number of dimensions
mx = size(x,1);                    % number of trial sites and their dimension
% normalize trial sites
x = (x - repmat(dmodel.Ssc(1,:),mx,1))./repmat(dmodel.Ssc(2,:),mx,1);
%  get distances to design sites
dx = zeros(mx*m,n);  kk = 1:m;
for  k = 1 : mx
    dx(kk,:) = repmat(x(k,:),size(dmodel.S,1),1) - dmodel.S;
    kk = kk + m;
end

f = feval(dmodel.regr, x);
r= feval(dmodel.corr, dmodel.theta, dx);
r = reshape(r, m, mx);
%----------------------------predictor
% Scaled predictor
sy = f * dmodel.beta + (dmodel.gamma * r)';
%  Predictor
y = dmodel.Ysc(1,:) + dmodel.Ysc(2,:) * sy;
%----------------------------mean square error
rt = dmodel.C \ r;
u = dmodel.G \ (dmodel.Ft.' * rt - f.');
mse = dmodel.sigma2 .* (1 + u.^2 - sum(rt.^2,1))';





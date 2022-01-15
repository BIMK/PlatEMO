function K = covLIN(hyp, x, z, i)

% Linear covariance function. The covariance function is parameterized as:
%
% k(x^p,x^q) = x^p'*x^q
%
% The are no hyperparameters:
%
% hyp = [ ]
%
% Note that there is no bias or scale term; use covConst to add these.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '0'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

% compute inner products
if dg                                                               % vector kxx
  K = sum(x.*x,2);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = x*x';
  else                                                   % cross covariances Kxz
    K = x*z';
  end
end

if nargin>3                                                        % derivatives
  error('Unknown hyperparameter')
end
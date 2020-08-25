function [e, g] = gp_eg(w, gp, x, y, varargin)
%GP_EG  Evaluate the energy function (un-normalized negative marginal
%       log posterior) and its gradient
%
%  Description
%    [E, G] = GP_EG(W, GP, X, Y, OPTIONS) takes a Gaussian process
%    structure GP together with a matrix X of input vectors and a
%    matrix Y of targets, and evaluates the energy function E and
%    its gradient G. Each row of X corresponds to one input vector
%    and each row of Y corresponds to one target vector.
%
%    The energy is minus log posterior cost function:
%        E = EDATA + EPRIOR 
%          = - log p(Y|X, th) - log p(th),
%    where th represents the parameters (lengthScale,
%    magnSigma2...), X is inputs and Y is observations (regression)
%    or latent values (non-Gaussian likelihood).
%
%    OPTIONS is optional parameter-value pair
%      z - optional observed quantity in triplet (x_i,y_i,z_i)
%          Some likelihoods may use this. For example, in case of
%          Poisson likelihood we have z_i=E_i, that is, expected
%          value for ith case.
%
%  See also
%    GP_E, GP_G
%
% Copyright (c) 2010 Aki Vehtari
  
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% Single function for some optimization routines, no need for mydeal...
e=gp_e(w, gp, x, y, varargin{:});
if nargout>1
  if isnan(e)
    g=NaN;
  else
    g=gp_g(w, gp, x, y, varargin{:});
  end
end

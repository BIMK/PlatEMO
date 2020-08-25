function gp = gp_unpak(gp, w, param)
%GP_UNPAK  Set GP parameters from vector to structure
%
%  Description
%    GP = GP_UNPAK(GP, W, PARAM) takes an Gaussian Process data
%    structure GP and a parameter vector W, and returns a Gaussian
%    Process structure identical to the input, except that the
%    parameters has been set to the ones in W. PARAM defines which
%    parameters are present in the W vector. If PARAM is not given
%    the function unpacks all parameters.
%
%    Each of the following strings in PARAM defines one group of
%    parameters to unpack:
%      covariance  - unpack parameters of covariance function
%      likelihood  - unpack parameters of likelihood
%      inducing    - unpack inducing inputs (in sparse approximations): 
%                    W = gp.X_u(:)
%
%    By combining the strings one can unpack more than one group of
%    parameters. For example:
%      covariance+inducing   - unpack covariance function parameters
%                              and inducing inputs
%      covariance+likelihood - unpack covariance function parameters
%                              and likelihood parameters
%
%    Inside each group (such as covariance functions) the
%    parameters to be unpacked is defined by the existence of a
%    prior structure. For example, if GP has two covariance
%    functions but only the first one has prior for its parameters
%    then only the parameters of the first one are unpacked. Thus,
%    also inducing inputs require prior if they are to be
%    optimized.
% 
%    GP_PAK and GP_UNPAK functions are used, e.g., when GP
%    parameters are optimized with GP_OPTIM or sampled with GP_MC. 
%    See GP_SET and option 'infer_params'.
%
%  See also
%    GP_PAK, GP_SET
%
% Copyright (c) 2007-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if isempty(w)
  return
end

if nargin < 3
  param = gp.infer_params;
end

if size(w,1) > 1
  error(' The vector to be packed has to be row vector! \n')
end

% Unpack the parameters of covariance functions
if ~isempty(strfind(param, 'covariance'))
  ncf = length(gp.cf);
  
  for i=1:ncf
    gpcf = gp.cf{i};
    [gpcf, w] = gpcf.fh.unpak(gpcf, w);
    gp.cf{i} = gpcf;
  end
  
end

% Unpack the parameters of likelihood function
if ~isempty(strfind(param, 'likelihood'))
  [gp.lik w] = gp.lik.fh.unpak(gp.lik, w);
  
  % Unpack the parameters of the second likelihood function (monotonicity)
  if isfield(gp, 'lik_mono')
    [gp.lik_mono w] = gp.lik_mono.fh.unpak(gp.lik_mono, w);
  end
end


% Unpack the inducing inputs
if ~isempty(strfind(param, 'inducing'))
  if isfield(gp,'p') && isfield(gp.p, 'X_u') && ~isempty(gp.p.X_u)
    if ~iscell(gp.p.X_u)
      % One prior for all inducing inputs
      lu = length(gp.X_u(:));
      gp.X_u = reshape(w(1:lu), size(gp.X_u'))';
      if lu < length(w)
        w = w(lu+1:end);
      end
      [gp.p.X_u, w]=gp.p.X_u.fh.unpak(gp.p.X_u, w);
    else
      % Own prior for each inducing input
      d=numel(gp.X_u(1,:));
      for i=1:size(gp.X_u)
        gp.X_u(i,:) = reshape(w(1:d), size(gp.X_u(1,:)));
        w(1:d)=[];
        [gp.p.X_u{i}, w] = gp.p.X_u{i}.fh.unpak(gp.p.X_u{i}, w);
      end
    end
  end
end

% Unpack the covariance function parameters
if ~isempty(strfind(param, 'mean'))
  nmf = length(gp.meanf);
  
  for i=1:nmf
    gpmf = gp.meanf{i};
    [gpmf, w] = gpmf.fh.unpak(gpmf, w);
    gp.meanf{i} = gpmf;
  end
  
end

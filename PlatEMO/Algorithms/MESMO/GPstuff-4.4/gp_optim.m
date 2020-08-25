function [gp, varargout] = gp_optim(gp, x, y, varargin)
%GP_OPTIM  Optimize paramaters of a Gaussian process 
%
%  Description
%    GP = GP_OPTIM(GP, X, Y, OPTIONS) optimises the parameters of a
%    GP structure given matrix X of training inputs and vector
%    Y of training targets.
%
%    [GP, OUTPUT1, OUTPUT2, ...] = GP_OPTIM(GP, X, Y, OPTIONS)
%    optionally returns outputs of the optimization function.
%
%    OPTIONS is optional parameter-value pair
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of
%               Poisson likelihood we have z_i=E_i, that is, expected
%               value for ith case.
%      optimf - function handle for an optimization function, which is
%               assumed to have similar input and output arguments
%               as usual fmin*-functions. Default is @fminscg.
%      opt    - options structure for the minimization function. 
%               Use optimset to set these options. By default options
%               'GradObj' is 'on', 'LargeScale' is 'off'.
%      loss   - 'e' to minimize the marginal posterior energy (default) or
%               'loo' to minimize the negative leave-one-out lpd
%               'kfcv' to minimize the negative k-fold-cv lpd
%               'waic' to minimize the WAIC loss
%               only 'e' and 'loo' with Gaussian likelihood have gradients
%      k      - number of folds in kfcv
%
%  See also
%    GP_SET, GP_E, GP_G, GP_EG, FMINSCG, FMINLBFGS, OPTIMSET, DEMO_REGRESSION*
%
% Copyright (c) 2010-2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.FunctionName = 'GP_OPTIM';
ip.addRequired('gp',@(x) isstruct(x) || isempty(x));
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('optimf', @fminscg, @(x) isa(x,'function_handle'))
ip.addParamValue('opt', [], @isstruct)
ip.addParamValue('loss', 'e', @(x) ismember(lower(x),{'e', 'loo', 'kfcv', 'waic' 'waic' 'waicv' 'waicg'}))
ip.addParamValue('k', 10, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
ip.parse(gp, x, y, varargin{:});
if isempty(gp)
  gp=gp_set();
end
if isempty(gp_pak(gp))
  % nothing to optimize
  return
end
z=ip.Results.z;
optimf=ip.Results.optimf;
opt=ip.Results.opt;
loss=ip.Results.loss;
k=ip.Results.k;


switch lower(loss)
  case 'e'
    fh_eg=@(ww) gp_eg(ww, gp, x, y, 'z', z);
    optdefault=struct('GradObj','on','LargeScale','off');
  case 'loo'
    fh_eg=@(ww) gp_looeg(ww, gp, x, y, 'z', z);
    if isfield(gp.lik.fh,'trcov') || isequal(gp.latent_method, 'EP')
      optdefault=struct('GradObj','on','LargeScale','off');
    else
      % Laplace-LOO does not have yet gradients
      optdefault=struct('Algorithm','interior-point');
      if ismember('optimf',ip.UsingDefaults)
        optimf=@fmincon;
      end
    end
  case 'kfcv'
    % kfcv does not have yet gradients
    fh_eg=@(ww) gp_kfcve(ww, gp, x, y, 'z', z, 'k', k);
    optdefault=struct('Algorithm','interior-point');
    if ismember('optimf',ip.UsingDefaults)
      optimf=@fmincon;
    end
  case {'waic' 'waicv'}
    % waic does not have yet gradients
    fh_eg=@(ww) -gp_waic(gp_unpak(gp,ww), x, y, 'z', z);
    optdefault=struct('Algorithm','interior-point');
    if ismember('optimf',ip.UsingDefaults)
      optimf=@fmincon;
    end
  case 'waicg'
    % waic does not have yet gradients
    fh_eg=@(ww) -gp_waic(gp_unpak(gp,ww), x, y, 'z', z, 'method', 'G');
    optdefault=struct('Algorithm','interior-point');
    if ismember('optimf',ip.UsingDefaults)
      optimf=@fmincon;
    end
end
opt=setOpt(optdefault,opt);
w=gp_pak(gp);
if isequal(lower(loss),'e') || (isequal(lower(loss),'loo')) && (isfield(gp.lik.fh,'trcov') || isequal(gp.latent_method, 'EP'))
  switch nargout
    case 6
      [w,fval,exitflag,output,grad,hessian] = optimf(fh_eg, w, opt);
      varargout={fval,exitflag,output,grad,hessian};
    case 5
      [w,fval,exitflag,output,grad] = optimf(fh_eg, w, opt);
      varargout={fval,exitflag,output,grad};
    case 4
      [w,fval,exitflag,output] = optimf(fh_eg, w, opt);
      varargout={fval,exitflag,output};
    case 3
      [w,fval,exitflag] = optimf(fh_eg, w, opt);
      varargout={fval,exitflag};
    case 2
      [w,fval] = optimf(fh_eg, w, opt);
      varargout={fval};
    case 1
      w = optimf(fh_eg, w, opt);
      varargout={};
  end
else
  lb=repmat(-8,size(w));
  ub=repmat(10,size(w));
  switch nargout
    case 6
      [w,fval,exitflag,output,grad,hessian] = optimf(fh_eg, w, [], [], [], [], lb, ub, [], opt);
      varargout={fval,exitflag,output,grad,hessian};
    case 5
      [w,fval,exitflag,output,grad] = optimf(fh_eg, w, [], [], [], [], lb, ub, [], opt);
      varargout={fval,exitflag,output,grad};
    case 4
      [w,fval,exitflag,output] = optimf(fh_eg, w, [], [], [], [], lb, ub, [], opt);
      varargout={fval,exitflag,output};
    case 3
      [w,fval,exitflag] = optimf(fh_eg, w, [], [], [], [], lb, ub, [], opt);
      varargout={fval,exitflag};
    case 2
      [w,fval] = optimf(fh_eg, w, [], [], [], [], lb, ub, [], opt);
      varargout={fval};
    case 1
      w = optimf(fh_eg, w, [], [], [], [], lb, ub, [], opt);
      varargout={};
  end
end
gp=gp_unpak(gp,w);
end

function opt=setOpt(optdefault, opt)
  % Set default options
  opttmp=optimset(optdefault,opt);
  
  % Set some additional options for @fminscg
  if isfield(opt,'lambda')
    opttmp.lambda=opt.lambda;
  end
  if isfield(opt,'lambdalim')
    opttmp.lambdalim=opt.lambdalim;
  end
  opt=opttmp;
end


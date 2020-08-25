function gp = gp_set(varargin)
%GP_SET  Create a Gaussian process model structure. 
%
%  Description
%    GP = GP_SET('PARAM1',VALUE1,'PARAM2',VALUE2,...)
%    creates a Gaussian process structure in which the named
%    parameters have the specified values. Any unspecified
%    parameters are set to default values. Either 'cf' or 
%    'meanf' parameter has to be specified. 
%  
%    GP = GP_SET(GP,'PARAM1',VALUE1,'PARAM2',VALUE2,...) 
%    modify a Gaussian process structure with the named
%    parameters altered with the specified values.
%
%    Parameters for Gaussian process
%      cf           - Single covariance structure or cell array of 
%                     covariance function structures created by
%                     gpcf_* functions. The default is [],
%                     except if neither cf or meanf is given then the
%                     the default is {gpcf_constant() gpcf_sexp()}. 
%      meanf        - Single mean function structure or cell array of 
%                     mean function function structures created by
%                     gpmf_* functions. The default is [].
%                     Mean functions work only with GP type 'FULL'
%      type         - Type of Gaussian process
%                      'FULL'   full GP (default)
%                      'FIC'    fully independent conditional sparse
%                               approximation
%                      'PIC'    partially independent conditional  
%                               sparse approximation
%                      'CS+FIC' compact support + FIC model sparse 
%                               approximation
%                      'DTC'    deterministic training conditional 
%                               sparse approximation
%                      'SOR'    subset of regressors sparse
%                               approximation
%                      'VAR'    variational sparse approximation
%      lik          - Likelihood structure created by one of the 
%                     likelihood functions lik_*. The default is
%                     created by lik_gaussian(). If likelihood is
%                     non-Gaussian, see latent_method below.
%      jitterSigma2 - Positive jitter to be added in the diagonal of 
%                     covariance matrix. The default is 0.
%      infer_params - String defining which parameters are inferred.
%                     The default is 'covariance+likelihood'.
%                      'covariance'     = infer parameters of the
%                                         covariance functions
%                      'likelihood'     = infer parameters of the likelihood
%                      'inducing'       = infer inducing inputs (in sparse
%                                         approximations): W = gp.X_u(:)    
%                       By combining the strings one can infer more than 
%                       one group of parameters. For example:
%                      'covariance+likelihood' = infer covariance function
%                                                and likelihood parameters
%                      'covariance+inducing' = infer covariance function
%                                              parameters and inducing 
%                                              inputs
%      comp_cf      - Option for when multiple latents are inferred
%                     (gpla2_e/g/pred). Sets specific covariance functions
%                     for specific latent processes. Must be given in cell
%                     format where ith element of the cell is vector
%                     defining covariance functions for ith latent process.
%      savememory   - Option for memory saving. Used in gradient
%                     calculations. The defaults is 'off'.
%
%    The additional fields when the likelihood is not Gaussian
%    (lik is not lik_gaussian or lik_gaussiansm) are:
%      latent_method - Method for marginalizing over latent
%                      values. Possible methods are 
%                      'Laplace' (default), 'EP' and 'MCMC'.
%      latent_opt    - Additional option structure for the chosen
%                      latent method. See default values for
%                      options below.
%    The options which can be set for each latent method are
%      MCMC:
%        method - Function handle to function which samples the
%                 latent values @esls (default), @scaled_mh or @scaled_hmc
%        f      - 1xn vector of latent values. The default is [].
%      Laplace:
%        optim_method - Method to find the posterior mode
%                      'newton' (default except for lik_t)
%                      'stabilized-newton', 'fminuc_large', or
%                      'lik_specific' (applicable and default for lik_t)
%        tol          - Iterations are stopped when change in the log 
%                       marginal likelihood is smaller than tol. The 
%                       default is 1e-6.
%      EP: 
%        maxiter      - Maximum number of EP iterations. The default is 20.
%        tol          - Iterations are stopped when all changes in the log 
%                       predictive densities and the log marginal likelihood 
%                       are smaller than tol. The default is 1e-4.
%        parallel     - Use parallel updating of site parameters: 
%                       'on' (default) or 'off'
%        init_prev    - Use parameter values from previous EP-iterations as
%                       initial parameter values: 'on' (default) or 'off'       
%        df           - Damping factor. Default is 0.8 for parallel-EP and 1.0 
%                       for sequential-EP.
%        optim_method - Method for evaluating EP. Default is 'basic-EP' for log
%                       concave likelihoods and 'robust-EP' for not log concave.                
%
%      for robust-EP
%        ninit        - Number of initial parallel iterations. Default is 10.
%        maxiter      - Maximum number of EP iterations. The default is 200.
%        df           - Initial damping factor. Default is 0.8.
%        eta          - Eta parameter for fractional EP. Default is 1.
%        eta2         - Secondary eta parameter for fractional EP, which is
%                       used only when the double-loop method is
%                       unable to proceed with eta. Default is 0.5.
%        max_ninner   - The maximum number of inner-loop iterations in the
%                       double-loop algorithm. The higher the
%                       value, the more robust the convergence
%                       properties are but with the cost of
%                       increased computational burden. Default is 4.
%        tolStop      - Tolerance level for stopping the iterations.
%                       Default is 1e-4.
%        tolGrad      - The minimum decrease gradient (g) in the search 
%                       direction, abs(g_new)<tolGrad*abs(g) This
%                       can be used to control the amount of step
%                       size adjustments Default is 0.9.
%        tolInner     - The inner loop energy tolerance. Smaller tolerance 
%                       increases the robustness, but increases the 
%                       computational cost. Default is 1e-3.
%        tolUpdate    - Tolerance level for ignoring updates. Default is 1e-6.
%        cavity_var_lim - Limit the for cavity variance Vc,
%                         proportional to the prior variance diag(K):
%                         Vc < Vc_lim*diag(K). Default is 2.
%        up_mode      - The search direction in the double-loop algorithm:
%                       'ep' or 'grad'. Default is 'ep'.
%                       Note that the implementation will
%                       automatically switch to gradients if the
%                       inner-loop minimization fails to to reduce
%                       the gradient within tolGrad
%        df_lim       - Limit for the step size. Default is 1.
%        display      - Control the amount of diagnostic verbosity.
%                       'off' displays nothing (default), 'final'
%                       display the hyperparameters and the final
%                       output, and 'iter' displays output at each
%                       iteration.
%        
%    The additional fields needed in sparse approximations are:
%      X_u          - Inducing inputs, no default, has to be set when
%                     FIC, PIC, PIC_BLOCK, VAR, DTC, or SOR is used.
%      Xu_prior     - Prior for inducing inputs. The default is prior_unif.
%
%    The additional field required by PIC sparse approximation is:
%      tr_index     - The blocks for the PIC model. The value has to
%                     be a cell array of the index vectors appointing
%                     the data points into blocks. For example, if x  
%                     is a matrix of data inputs then x(tr_index{i},:) 
%                     are the inputs belonging to the i'th block.
%
%    The additional fields needed with derivative observations
%      derivobs     - Tells whether derivative observations are
%                     used: 'on' or 'off' (default).
%
%  See also
%    GPCF_*, LIK_*, PRIOR_*, GP_PAK, GP_UNPAK, GP_E, GP_G, GP_EG,
%    GP_PRED, GP_MC, GP_IA, ...
%
%  References:
%    Quinonero-Candela, J. and Rasmussen, C. E. (2005). A unifying
%    view of sparse approximate Gaussian process regression. 
%    Journal of Machine Learning Research, 6(3):1939-1959.
%
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.
%
%    Snelson, E. and Ghahramani, Z. (2006). Sparse Gaussian process
%    using pseudo-inputs. In Weiss, Y., Sch�lkopf, B., and Platt,
%    J. (eds) Advances in Neural Information Processing Systems 18,
%    pp. 1257-1264.
%
%    Titsias, M. K. (2009). Variational Model Selection for Sparse
%    Gaussian Process Regression. Technical Report, University of
%    Manchester.
%
%    Vanhatalo, J. and Vehtari, A. (2008). Modelling local and
%    global phenomena with sparse Gaussian processes. Proceedings
%    of the 24th Conference on Uncertainty in Artificial
%    Intelligence.
%
%    Sarkka, S., Solin, A., Hartikainen, J. (2013). 
%    Spatiotemporal learning via infinite-dimensional Bayesian 
%    filtering and smoothing. IEEE Signal Processing Magazine, 
%    30(4):51-61.
%
% Copyright (c) 2006-2010 Jarno Vanhatalo
% Copyright (c) 2010-2011 Aki Vehtari
% Copyright (c) 2014 Arno Solin and Jukka Koskenranta
  
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GP_SET';
  ip.addOptional('gp', [], @isstruct);
  ip.addParamValue('cf',[], @(x) isempty(x) || isstruct(x) || iscell(x));
  ip.addParamValue('meanf',[], @(x) isempty(x) || isstruct(x) || iscell(x));
  ip.addParamValue('type','FULL', ...
                   @(x) ismember(x,{'FULL' 'FIC' 'PIC' 'PIC_BLOCK' 'VAR' ...
                      'DTC' 'SOR' 'CS+FIC','KALMAN'}));
  ip.addParamValue('lik',lik_gaussian(), @(x) isstruct(x));
  ip.addParamValue('jitterSigma2',0, @(x) isscalar(x) && x>=0);
  ip.addParamValue('infer_params','covariance+likelihood', @(x) ischar(x));
  ip.addParamValue('latent_method','Laplace', @(x) ischar(x) || iscell(x));
  ip.addParamValue('latent_opt',struct(), @isstruct);
  ip.addParamValue('X_u',[],  @(x) isreal(x) && all(isfinite(x(:))));
  ip.addParamValue('Xu_prior',prior_unif,  @(x) isstruct(x) || isempty(x) || ...
                   iscell(x));
  ip.addParamValue('tr_index', [], @(x) ~isempty(x) || iscell(x))    
  ip.addParamValue('comp_cf', [], @(x) iscell(x))    
  ip.addParamValue('derivobs','off', @(x) islogical(x) || isscalar(x) || ...
                   (ischar(x) && ismember(x,{'on' 'off'})));
  ip.addParamValue('savememory','off', @(x) islogical(x) || isscalar(x) || ...
                   (ischar(x) && ismember(x,{'on' 'off'})));
%   ip.addParamValue('optim_method', [], @(x) isreal(x) && (x==1 || x==2) &&  ...
%                     isfinite(x))
  ip.parse(varargin{:});
  gp=ip.Results.gp;


  if isempty(gp)
    % Initialize a Gaussian process
    init=true;
  else
    % Modify a Gaussian process
    if ~isfield(gp,'cf')
      error('First argument does not seem to be a Gaussian process structure')
    end
    init=false;
  end

  % FULL or sparse
  if init || ~ismember('type',ip.UsingDefaults)
    gp.type=ip.Results.type;
  end
  % Likelihood
  if init || ~ismember('lik',ip.UsingDefaults)
    gp.lik = ip.Results.lik;
  end
  % Covariance function(s)
  if init || ~ismember('cf',ip.UsingDefaults)
    gp.cf=ip.Results.cf;
    if isstruct(gp.cf)
      % store single structure in a cell array, too
      gp.cf={gp.cf};
    end
  end
  % Mean function(s)
  if init || ~ismember('meanf',ip.UsingDefaults)
    if ~isempty(ip.Results.meanf)
      if ~isequal(gp.type,'FULL')
        error('Mean functions ''meanf'' can be used only with GP type ''FULL''');
      end
      gp.meanf=ip.Results.meanf;
      if isstruct(gp.meanf)
        % store single structure in a cell array, too
        gp.meanf={gp.meanf};
      end
    end
  end
  if isempty(gp.cf) && (~isfield(gp,'meanf') || isempty(gp.meanf))
    gp.cf={gpcf_constant() gpcf_sexp()};
  end
  % Inference for which parameters 
  if init || ~ismember('infer_params',ip.UsingDefaults)
    gp.infer_params=ip.Results.infer_params;
  end
  % Jitter
  if init || ~ismember('jitterSigma2',ip.UsingDefaults)
    gp.jitterSigma2=ip.Results.jitterSigma2;
  end
  % Gradient observation
  if init || ~ismember('derivobs',ip.UsingDefaults) || ~isfield(gp,'derivobs')
    derivobs=ip.Results.derivobs;
    if ~ischar(derivobs)
      if derivobs
        derivobs='on';
      else
        derivobs='off';
      end
    end
    switch derivobs
      case 'on'
        gp.derivobs=true;
      case 'off'
        if isfield(gp,'derivobs')
          gp=rmfield(gp,'derivobs');
        end
    end
  end
  % Covariance functions for multiple latent processes
  if ~ismember('comp_cf',ip.UsingDefaults)
    gp.comp_cf=ip.Results.comp_cf;
  end
  
  if init || ~ismember('savememory',ip.UsingDefaults) || ~isfield(gp,'savememory')
    savememory=ip.Results.savememory;
    if ~ischar(savememory)
      if savememory
        savememory='on';
      else
        savememory='off';
      end
    end
    switch savememory
      case 'on'
        gp.savememory=true;
      case 'off'
        if isfield(gp,'savememory')
          gp=rmfield(gp,'savememory');
        end
    end
  end

  % Inducing inputs
  if ismember(gp.type,{'FIC' 'CS+FIC' 'DTC' 'VAR' 'SOR' 'PIC' 'PIC_BLOCK'})
    if init || ~ismember('X_u',ip.UsingDefaults)
      gp.X_u = ip.Results.X_u;
      gp.nind = size(gp.X_u,1);
    end
    if init || ~ismember('Xu_prior',ip.UsingDefaults) || ~isfield(gp,'p')
      gp.p.X_u = ip.Results.Xu_prior;
    end
    if ismember(gp.type, {'PIC' 'PIC_BLOCK'})
      % + PIC block indexes
      if init || ~ismember('tr_index',ip.UsingDefaults)
        gp.tr_index = ip.Results.tr_index;
      end
    end
  end
  if ismember(gp.type,{'FIC' 'PIC' 'PIC_BLOCK' 'VAR' 'DTC' 'SOR'}) ...
      && isempty(gp.X_u)
    error(sprintf('Need to set X_u when using %s',gp.type))
  end
  if ismember(gp.type,{'CS+FIC'})
    % check that we have both cs and non-cs type covariance functions
    ncf=numel(gp.cf);
    iscs=0;isnoncs=0;
    for i = 1:ncf
      if isfield(gp.cf{i},'cs')
        iscs=1;
      else
        isnoncs=1;
      end
    end
    if ~(iscs && isnoncs)
      error('With CS+FIC need to define at least one cs and one non-cs covariance function')
    end
  end
  if ismember(gp.type,{'KALMAN'})
      
    % Check likelihood function
    if ~strcmpi(gp.lik.type,'Gaussian')
      error('The ''KALMAN'' option only supports Gaussian likelihoods.')
    end
    
    % Check covariance functions
    for j=1:length(gp.cf)
      if ~isfield(gp.cf{j}.fh,'cf2ss'),
         error('State space conversion not implemented for ''%s''.', ...
             gp.cf{j}.type) 
      end
    end
    
  end
  % Latent method
  if isfield(gp.lik.fh,'trcov') && ~isfield(gp, 'lik_mono')
    % Gaussian likelihood
    if ~ismember('latent_method',ip.UsingDefaults)
      error('No latent method needed with a Gaussian likelihood')
    end
    if isfield(gp,'latent_method')
      gp=rmfield(gp,'latent_method')
    end
    gp.fh.e=@gp_e;
    gp.fh.g=@gp_g;
    gp.fh.pred=@gp_pred;
    gp.fh.jpred=@gp_jpred;
  else
    if init || ~ismember('latent_method',ip.UsingDefaults) || ~isfield(gp,'latent_method')
      latent_method=ip.Results.latent_method;
      switch latent_method
        case 'MCMC'
          % Remove traces of other latent methods
          if isfield(gp,'latent_opt'); gp=rmfield(gp,'latent_opt'); end
          if isfield(gp,'fh') && isfield(gp.fh,'ne')
            gp.fh=rmfield(gp.fh,{'ne' 'e' 'g' 'pred' 'jpred' 'looe' 'loog'}); 
          end
          % Set latent method
          gp.latent_method=latent_method;
          gp.fh.pred=@gpmc_pred;
          gp.fh.jpred=@gpmc_jpred;
        case 'EP'
          % Remove traces of other latent methods
          if isfield(gp,'latent_method') && ~isequal(latent_method,gp.latent_method) && isfield(gp,'latent_opt')
            gp=rmfield(gp,'latent_opt');
          end
          if isfield(gp,'latentValues'); gp=rmfield(gp,'latentValues'); end
          % Set latent method
          gp.latent_method=latent_method;
          % following sets gp.fh.e = @ep_algorithm;
          gp = gpep_e('init', gp);
        case 'Laplace'
          % Remove traces of other latent methods
          if isfield(gp,'latent_method') && ~isequal(latent_method,gp.latent_method) && isfield(gp,'latent_opt')
            gp=rmfield(gp,'latent_opt');
          end
          if isfield(gp,'latentValues'); gp=rmfield(gp,'latentValues'); end
          % Set latent method
          gp.latent_method=latent_method;
          % following sets gp.fh.e = @laplace_algorithm;
          gp = gpla_e('init', gp);
        case 'NA'
          % no latent method set
          if isfield(gp,'latent_method'); gp=rmfield(gp,'latent_method'); end
          if isfield(gp,'latent_opt'); gp=rmfield(gp,'latent_opt'); end
        otherwise
          error('Unknown type of latent_method!')
      end % switch latent_method
    end % if init || ~ismember('latent_method',ip.UsingDefaults)
    if init || ~ismember('latent_opt',ip.UsingDefaults) || ~isfield(gp,'latent_opt')
      latent_opt=ip.Results.latent_opt;
      switch gp.latent_method
        case 'MCMC'
          % Handle latent_opt
          ipmc=inputParser;
          ipmc.FunctionName = 'GP_SET - latent method MCMC options';
          ipmc.addParamValue('method',@esls, @(x) isa(x,'function_handle'));
          ipmc.addParamValue('f',[],  @(x) isreal(x) && all(isfinite(x(:))));
          ipmc.parse(latent_opt);
          if init || ~ismember('method',ipmc.UsingDefaults) || ~isfield(gp.fh,'mc')
            gp.fh.mc = ipmc.Results.method;
          end
          if init || ~ismember('f',ipmc.UsingDefaults) || ~isfield(gp,'latentValues')
            gp.latentValues = ipmc.Results.f;
          end
        case 'EP'
          % Handle latent_opt
          ipep=inputParser;
          ipep.FunctionName = 'GP_SET - latent method EP options';
          ipep.addParamValue('optim_method',[], @(x) ischar(x));
          ipep.addParamValue('maxiter',20, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0);
          ipep.addParamValue('display', 'off', @(x) ischar(x) && ismember(x,{'off', 'final', 'iter'}))
          ipep.addParamValue('parallel','on', @(x) ischar(x) && ismember(x,{'off', 'on'}));    % default on
          ipep.addParamValue('init_prev', 'on', @(x) ischar(x) && ismember(x,{'off', 'on'}));    % default on
          % Following option is only for basic-EP
          % all changes in the log predictive densities and the log marginal
          % likelihood are smaller than tol.
          ipep.addParamValue('tol',1e-4, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0);
          % Following options are only for robust-EP
          % max number of initial parallel iterations
          ipep.addParamValue('ninit', 10, @(x) isreal(x) && (x==1 || rem(1,x)==1) && isfinite(x))
          % max number of inner loop iterations in the double-loop algorithm
          ipep.addParamValue('max_ninner', 4, @(x) isreal(x) && (x==1 || rem(1,x)==1) && isfinite(x))
          % converge tolerance with respect to the maximum change in E[f] and Var[f]
          ipep.addParamValue('tolStop', 1e-4, @(x) isreal(x) && isfinite(x))
          % tolerance for the EP site updates
          ipep.addParamValue('tolUpdate', 1e-6, @(x) isreal(x) && isfinite(x))
          % inner loop energy tolerance
          ipep.addParamValue('tolInner', 1e-3, @(x) isreal(x) && isfinite(x))
          % minimum gradient (g) decrease in the search direction, abs(g_new)<tolGrad*abs(g)
          ipep.addParamValue('tolGrad', 0.9, @(x) isreal(x) && isfinite(x))
          % limit for the cavity variance Vc, Vc < Vc_lim*diag(K)
          ipep.addParamValue('cavity_var_lim', 2, @(x) isreal(x) && isfinite(x))
          % the intial damping factor
          ipep.addParamValue('df', 0.8, @(x) isreal(x) && isfinite(x))
          % the initial fraction parameter
          ipep.addParamValue('eta', 1, @(x) isreal(x) && isfinite(x))
          % the secondary fraction parameter          
          ipep.addParamValue('eta2', .5, @(x) isreal(x) && isfinite(x))
          % update mode in double-loop iterations
          ipep.addParamValue('up_mode', 'ep', @(x) ischar(x) && ismember(x,{'ep' 'grad'}))
          % step size limit (1 suitable for ep updates)
          ipep.addParamValue('df_lim', 1, @(x) isreal(x) && isfinite(x))
          % whether to do one inner loop-iteration per site ('on') or until
          % convergence ('off') in nester ep
          ipep.addParamValue('incremental', 'on', @(x) ischar(x) && ismember(x,{'off', 'on'}));    % default on                    
          ipep.parse(latent_opt);
          optim_method = ipep.Results.optim_method;
          if ~isempty(optim_method)
            gp.latent_opt.optim_method=optim_method;
          else
            % If likelihood is not log-concave (exists functions siteDeriv2
            % & tiltedMoments2) use robust-EP by default, else basic EP.
            if isfield(gp.lik.fh, 'siteDeriv2')
              gp.latent_opt.optim_method='robust-EP';
            else
              gp.latent_opt.optim_method='basic-EP';
            end
          end
          if init || ~ismember('maxiter',ipep.UsingDefaults) || ~isfield(gp,'maxiter')
            if strcmp(gp.latent_opt.optim_method,'robust-EP') && ismember('maxiter',ipep.UsingDefaults)
              gp.latent_opt.maxiter = 200;
            else
              gp.latent_opt.maxiter = ipep.Results.maxiter;
            end
          end
          if init || ~ismember('tol',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'tol')
            gp.latent_opt.tol = ipep.Results.tol;
          end
          if init || ~ismember('parallel',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'parallel')
            gp.latent_opt.parallel = ipep.Results.parallel;
          end
          if init || ~ismember('init_prev',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'init_prev')
            gp.latent_opt.init_prev = ipep.Results.init_prev;
          end
          if init || ~ismember('df',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'df')
            if strcmp(gp.latent_opt.parallel,'off') && ismember('df',ipep.UsingDefaults)
              gp.latent_opt.df = 1;
            else
              gp.latent_opt.df = ipep.Results.df;
            end
          end
          if strcmp(gp.latent_opt.optim_method, 'robust-EP')
            if init || ~ismember('ninit',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'ninit')
              gp.latent_opt.ninit = ipep.Results.ninit;
            end
            if init || ~ismember('eta',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'eta')
              gp.latent_opt.eta = ipep.Results.eta;
            end
            if init || ~ismember('eta2',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'eta2')
              gp.latent_opt.eta2 = ipep.Results.eta2;
            end
            if init || ~ismember('max_ninner',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'max_ninner')
              gp.latent_opt.max_ninner = ipep.Results.max_ninner;
            end
            if init || ~ismember('tolStop',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'tolStop')
              gp.latent_opt.tolStop = ipep.Results.tolStop;
            end
            if init || ~ismember('tolGrad',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'tolGrad')
              gp.latent_opt.tolGrad = ipep.Results.tolGrad;
            end
            if init || ~ismember('tolInner',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'tolInner')
              gp.latent_opt.tolInner = ipep.Results.tolInner;
            end
            if init || ~ismember('tolUpdate',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'tolUpdate')
              gp.latent_opt.tolUpdate = ipep.Results.tolUpdate;
            end
            if init || ~ismember('cavity_var_lim',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'cavity_var_lim')
              gp.latent_opt.cavity_var_lim = ipep.Results.cavity_var_lim;
            end
            if init || ~ismember('up_mode',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'up_mode')
              gp.latent_opt.up_mode = ipep.Results.up_mode;
            end
            if init || ~ismember('df_lim',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'df_lim')
              gp.latent_opt.df_lim = ipep.Results.df_lim;
            end
            if init || ~ismember('display',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'display')
              gp.latent_opt.display = ipep.Results.display;
            end
          end
          if isfield(gp.lik, 'nondiagW')
            if init || ~ismember('display',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'display')
              gp.latent_opt.display = ipep.Results.display;
            end
            if init || ~ismember('incremental',ipep.UsingDefaults) || ~isfield(gp.latent_opt,'incremental')
              gp.latent_opt.display = ipep.Results.display;
            end
          end
        case 'Laplace'
          % Handle latent_opt
          ipla=inputParser;
          ipla.FunctionName = 'GP_SET - latent method Laplace options';
          ipla.addParamValue('optim_method',[], @(x) ischar(x));
          ipla.addParamValue('tol',1e-4, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0);
          ipla.addParamValue('maxiter', 40, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0);
          ipla.parse(latent_opt);
          optim_method=ipla.Results.optim_method;
          if ~isempty(optim_method)
            gp.latent_opt.optim_method=optim_method;
          else
            if isfield(gp.lik.fh, 'optimizef')
              % slower than newton but more robust
              gp.latent_opt.optim_method='lik_specific'; 
            else
              gp.latent_opt.optim_method='newton';
            end
          end
          if init || ~ismember('tol',ipla.UsingDefaults) || ~isfield(gp.latent_opt,'tol')
            gp.latent_opt.tol = ipla.Results.tol;
          end
          if init || ~ismember('maxiter',ipla.UsingDefaults) || ~isfield(gp,'maxiter')
            gp.latent_opt.maxiter = ipla.Results.maxiter;
          end
        otherwise
          error('Unknown type of latent_method!')
      end % switch latent_method
    end % if init || ~ismember('latent_method',ip.UsingDefaults)
  end
  
end

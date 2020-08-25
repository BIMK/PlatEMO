function [w, s, h] = gp_pak(gp, param)
%GP_PAK  Combine GP parameters into one vector
%
%  Description
%    W = GP_PAK(GP, PARAM) takes a Gaussian Process structure
%    GP and string PARAM defining, which parameters are packed and
%    combines the parameters into a single row vector W. If PARAM
%    is not given the function packs all parameters.
%
%    Each of the following strings in PARAM defines one group of
%    parameters to pack:
%      covariance  - pack parameters of covariance function
%      likelihood  - pack parameters of likelihood
%      inducing    - pack inducing inputs (in sparse approximations): 
%                    W = gp.X_u(:)
%
%    By combining the strings one can pack more than one group of
%    parameters. For example:
%      covariance+inducing  - pack covariance function parameters
%                             and inducing inputs
%      covariance+likelih   - pack covariance function parameters
%                             and likelihood parameters
%
%    Inside each group (such as covariance functions) the
%    parameters to be packed is defined by the existence of a prior
%    structure. For example, if GP has two covariance functions but
%    only the first one has prior for its parameters then only the
%    parameters of the first one are packed. Thus, also inducing
%    inputs require prior if they are to be optimized.
%
%    GP_PAK and GP_UNPAK functions are used, e.g., when GP
%    parameters are optimized with GP_OPTIM or sampled with GP_MC.
%    See GP_SET and option 'infer_params'.
%
%    [W, WS] = GP_PAK(GP, PARAM) returns also cell array of string
%    labels for the weight vector elements, which makes diagnostics
%    easier.
%
%    [W, WS, H] = GP_PAK(GP, PARAM) returns also hierarchy level H of 
%    different parameters (0 for likelihood parameters and 1 for covariance
%    function parameters)
%
%  See also
%    GP_UNPAK, GP_SET
%
% Copyright (c) 2007-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  w = []; s = {}; h = [];

  if isfield(gp,'etr') && length(gp.etr) > 1
    if strcmp(gp.type, 'PIC_BLOCK') || strcmp(gp.type, 'PIC')
      ind = gp.tr_index;           % block indeces for training points
      gp = rmfield(gp,'tr_index');
    end
    ns = length(gp.etr);
    for i1 = 1:ns
      Gp = take_nth(gp,i1);
      [w(i1,:), s, h] = gp_pak(Gp);
    end
  else
    
    if nargin < 2
      param = gp.infer_params;
    end
    
    % Pack the parameters of covariance functions
    if ~isempty(strfind(param, 'covariance'))
      ncf = length(gp.cf);
      
      for i=1:ncf
        gpcf = gp.cf{i};
        [wi, si, hi] = gpcf.fh.pak(gpcf);
        w = [w wi];
        s = [s; si];
        h = [h hi];
      end
    end
    
    % Pack the parameters of likelihood function
    if ~isempty(strfind(param, 'likelihood'))
      [wi, si, hi] = gp.lik.fh.pak(gp.lik);
      w = [w wi];
      s = [s; si];
      h = [h hi];
    end
    
    % Pack the parameters of the second likelihood function (monotonic)
    if ~isempty(strfind(param, 'likelihood')) && isfield(gp, 'lik_mono')
      [wi, si, hi] = gp.lik_mono.fh.pak(gp.lik_mono);
      w = [w wi];
      s = [s; si];
      h = [h hi];
    end
    
    % Pack the inducing inputs
    if ~isempty(strfind(param, 'inducing'))
      if isfield(gp,'p') && isfield(gp.p, 'X_u') && ~isempty(gp.p.X_u)
        if ~iscell(gp.p.X_u)
          % One prior for all inducing inputs
          w = [w reshape(gp.X_u', numel(gp.X_u),1)'];
%           w = [w gp.X_u'];
          s = [s; sprintf('inducing x %d',numel(gp.X_u))];
          h = [h ones(1,numel(gp.X_u))];
          [wi,si,hi]=gp.p.X_u.fh.pak(gp.p.X_u);
          w = [w wi];
          s = [s; si];
          h = [h 1+hi];
        else
          % Own prior for each inducing input
          for i=1:size(gp.X_u,1)
            w = [w gp.X_u(i,:)];
            s = [s; sprintf('inducing x %d',numel(gp.X_u(i,:)))];
            h = [h ones(1,size(gp.X_u(i,:),2))];
            [wi,si,hi]=gp.p.X_u{i}.fh.pak(gp.p.X_u{i});
            si=strcat(repmat('prior-', size(si,1),1),si);
            w = [w wi];
            s = [s; si];
            h = [h 1+hi];            
          end
        end
      end
    end
    
    % Pack the prior weights and variances of mean functions
    if ~isempty(strfind(param, 'mean'))
      mf = length(gp.meanf);
      
      for i=1:mf
        gpmf = gp.meanf{i};
        [wi, si, hi] = gpmf.fh.pak(gpmf);
        w = [w wi];
        s = [s; si];
        h = [h, hi];
      end
    end
    
  end

end
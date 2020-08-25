function p = prior_gamma(varargin)
%PRIOR_GAMMA  Gamma prior structure     
%
%  Description
%    P = PRIOR_GAMMA('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) 
%    creates Gamma prior structure in which the named parameters
%    have the specified values. Any unspecified parameters are set
%    to default values.
%
%    P = PRIOR_GAMMA(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a prior structure with the named parameters altered
%    with the specified values.
%  
%    The parameterization is as in Gelman, Carlin, Stern, Dunson, Vehtari,
%    and Rubin (2013). Bayesian Data Analysis, third edition.
%
%    Parameters for Gamma prior [default]
%      sh       - shape (alpha) [4]
%      is       - inverse scale (beta) [1]
%      sh_prior - prior for sh [prior_fixed]
%      is_prior - prior for is [prior_fixed]
%
%  See also
%    PRIOR_*
%
% Copyright (c) 2000-2001,2010 Aki Vehtari
% Copyright (c) 2010 Jaakko Riihimäki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_GAMMA';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('sh',4, @(x) isscalar(x) && x>0);
  ip.addParamValue('sh_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('is',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('is_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Gamma';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Gamma')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('sh',ip.UsingDefaults)
    p.sh = ip.Results.sh;
  end
  if init || ~ismember('is',ip.UsingDefaults)
    p.is = ip.Results.is;
  end
  % Initialize prior structure
  if init
    p.p=[];
  end
  if init || ~ismember('sh_prior',ip.UsingDefaults)
    p.p.sh=ip.Results.sh_prior;
  end
  if init || ~ismember('is_prior',ip.UsingDefaults)
    p.p.is=ip.Results.is_prior;
  end

  if init
    % set functions
    p.fh.pak = @prior_gamma_pak;
    p.fh.unpak = @prior_gamma_unpak;
    p.fh.lp = @prior_gamma_lp;
    p.fh.lpg = @prior_gamma_lpg;
    p.fh.recappend = @prior_gamma_recappend;
  end

end

function [w, s, h] = prior_gamma_pak(p)
  
  w=[];
  s={};
  h=[];
  if ~isempty(p.p.sh)
    w = log(p.sh);
    s=[s; 'log(Gamma.sh)'];
    h = 1;
  end
  if ~isempty(p.p.is)
    w = [w log(p.is)];
    s=[s; 'log(Gamma.is)'];
    h = [h 1];
  end
end

function [p, w] = prior_gamma_unpak(p, w)

  if ~isempty(p.p.sh)
    i1=1;
    p.sh = exp(w(i1));
    w = w(i1+1:end);
  end
  if ~isempty(p.p.is)
    i1=1;
    p.is = exp(w(i1));
    w = w(i1+1:end);
  end
end

function lp = prior_gamma_lp(x, p)
  
  lp = sum(-p.is.*x + (p.sh-1).*log(x) +p.sh.*log(p.is)  -gammaln(p.sh));
  
  if ~isempty(p.p.sh)
    lp = lp + p.p.sh.fh.lp(p.sh, p.p.sh) + log(p.sh);
  end
  if ~isempty(p.p.is)
    lp = lp + p.p.is.fh.lp(p.is, p.p.is) + log(p.is);
  end
end

function lpg = prior_gamma_lpg(x, p)
  
  lpg = (p.sh-1)./x - p.is;
  
  if ~isempty(p.p.sh)
    lpgsh = (sum(-digamma1(p.sh) + log(p.is) + log(x)) + p.p.sh.fh.lpg(p.sh, p.p.sh)).*p.sh + 1;
    lpg = [lpg lpgsh];
  end
  if ~isempty(p.p.is)
    lpgis = (sum(p.sh./p.is+x) + p.p.is.fh.lpg(p.is, p.p.is)).*p.is + 1;
    lpg = [lpg lpgis];
  end
  
end

function rec = prior_gamma_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
  if ~isempty(p.p.sh)
    rec.sh(ri,:) = p.sh;
  end
  if ~isempty(p.p.is)
    rec.is(ri,:) = p.is;
  end
end    

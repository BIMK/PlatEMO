function p = prior_logunif(varargin)
%PRIOR_LOGUNIF  Uniform prior structure for the logarithm of the parameter
%       
%  Description
%    P = PRIOR_LOGUNIF creates uniform prior structure for the
%    logarithm of the parameter.
%    
%  See also
%    PRIOR_*
%
% Copyright (c) 2009 Jarno Vanhatalo
% Copyright (c) 2010 Jaakko Riihimäki
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_LOGUNIFORM';
  ip.addOptional('p', [], @isstruct);
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Log-Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Log-Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_logunif_pak;
    p.fh.unpak = @prior_logunif_unpak;
    p.fh.lp = @prior_logunif_lp;
    p.fh.lpg = @prior_logunif_lpg;
    p.fh.recappend = @prior_logunif_recappend;
  end
end

function [w, s,h] = prior_logunif_pak(p)
  w=[];
  s={};
  h=[];
end

function [p, w] = prior_logunif_unpak(p, w)
  w = w;
  p = p;
end

function lp = prior_logunif_lp(x, p)
  lJ = -log(x);    % log(1/x) log(|J|) of transformation
  lp = sum(lJ);
end

function lpg = prior_logunif_lpg(x, p)
  lJg = -1./x;     % gradient of log(|J|) of transformation
  lpg = lJg;
end

function rec = prior_logunif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
end

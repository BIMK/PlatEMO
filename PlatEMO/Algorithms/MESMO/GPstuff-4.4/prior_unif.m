function p = prior_unif(varargin)
%PRIOR_UNIF  Uniform prior structure     
%       
%  Description
%    P = PRIOR_UNIF creates uniform prior structure.
%    
%  See also
%    PRIOR_*
%
% Copyright (c) 2009 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_UNIFORM';
  ip.addOptional('p', [], @isstruct);
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_unif_pak;
    p.fh.unpak = @prior_unif_unpak;
    p.fh.lp = @prior_unif_lp;
    p.fh.lpg = @prior_unif_lpg;
    p.fh.recappend = @prior_unif_recappend;
  end
  
end

function [w, s, h] = prior_unif_pak(p, w)
  w=[];
  s={};
  h=[];
end

function [p, w] = prior_unif_unpak(p, w)
  w = w;
  p = p;
end

function lp = prior_unif_lp(x, p)
  lp = 0;
end

function lpg = prior_unif_lpg(x, p)
  lpg = zeros(size(x));
end

function rec = prior_unif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
end


function C = trcov(gpcf, x)
% TRCOV     Evaluate training covariance matrix for covariance function
%           This is a mex-function that is called from gpcf_*_trcov
%           functions.
%
%         Description    
%         K = TRCOV(GP, TX) takes in Gaussian process GP and matrix TX
%         that contains training input vectors to GP. Returns
%         noiseless covariance matrix K. Every element ij of K
%         contains covariance between inputs i and j in TX.
%
%         [K, C] = GP_TRCOV(GP, TX) returns also the noisy
%         covariance matrix C.
%
%         See also
%         GP_TRCOV, GPCF_SEXP -> GPCF_SEXP_TRCOV

% Copyright (c) 2008-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

    C = NaN; % If mex file is not compiled return NaN, in which case the
             % matrix is evaluated in gpcf_*_trcov

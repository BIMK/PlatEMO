function [srgtSRGT srgtSTT] = srgtsKRGFit(srgtOPT)
%Function srgtsKRGFit fits the specified kriging model using the DACE
%toolbox of Lophaven et al. (2002).
%
%    [srgtSRGT srgtSTT] = srgtsKRGFit(srgtOPT)
%
%srgtSRGT is the surrogate structure that contains the following fields:
%* KRG_DACEModel: DACE model, a struct with the elements:
%    * regr   : function handle to the regression model.
%    * corr   : function handle to the correlation function.
%    * theta  : correlation function parameters.
%    * beta   : generalized least squares estimate.
%    * gamma  : correlation factors.
%    * sigma2 : maximum likelihood estimate of the process variance.
%    * S      : scaled design sites.
%    * Ssc    : scaling factors for design arguments.
%    * Ysc    : scaling factors for design ordinates.
%    * C      : Cholesky factor of correlation matrix.
%    * Ft     : Decorrelated regression matrix.
%    * G      : From QR factorization: Ft = Q*G'.
%
%srgtSTT is the state structure that contains the following fields:
%* FIT_Fn       : function handle of the fitting function.
%* FIT_FnVal    : value of the loss function (after fitting).
%* KRG_DACEPerf : struct with DACE performance information:
%    * nv     : Number of evaluations of objective function.
%    * perf   : (q+2)*nv array, where q is the number of elements in theta,
%               and the columns hold current values of [theta;  psi(theta);
%               type]. type = 1, 2 or 3, indicate 'start', 'explore' or
%               'move.' A negative value for type indicates an uphill step.
%
%Example:
%     % basic information about the problem
%     myFN = @cos;  % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     npoints = 5;
%     X = linspace(designspace(1), designspace(2), npoints)';
%
%     % evaluate analysis function at X points
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     options = srgtsKRGSetOptions(X, Y);
%
%     [surrogate state] = srgtsKRGFit(options)
%
%     surrogate =
%
%     KRG_DACEModel: [1x1 struct]
%
%     state =
%
%     KRG_DACEPerf: [1x1 struct]
%
%REFERENCES:
%
%Lophaven SN, Nielsen HB, and Søndergaard J, DACE - A MATLAB Kriging
%Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%Modelling, Technical University of Denmark, 2002.
%Available at: http://www2.imm.dtu.dk/~hbn/dace/.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch func2str(srgtOPT.FIT_Fn)
    case 'dace_fit'
        srgtSTT = srgtsFitCreateState(srgtOPT);
        if isempty(srgtOPT.KRG_LowerBound) % no optimization for theta
            [srgtSRGT.KRG_DACEModel, srgtSTT.KRG_DACEPerf, srgtSTT.FIT_FnVal] = dace_fit(...
                srgtOPT.P, ...
                srgtOPT.T, ...
                srgtOPT.KRG_RegressionModel, ...
                srgtOPT.KRG_CorrelationModel, ...
                srgtOPT.KRG_Theta0);
        else
            [srgtSRGT.KRG_DACEModel, srgtSTT.KRG_DACEPerf, srgtSTT.FIT_FnVal] = dace_fit(...
                srgtOPT.P, ...
                srgtOPT.T, ...
                srgtOPT.KRG_RegressionModel, ...
                srgtOPT.KRG_CorrelationModel, ...
                srgtOPT.KRG_Theta0, ...
                srgtOPT.KRG_LowerBound, ...
                srgtOPT.KRG_UpperBound);
        end
        
    case 'srgtsXVFit'
        [srgtSRGT srgtSTT] = srgtsXVFit(srgtOPT);
        
end

return

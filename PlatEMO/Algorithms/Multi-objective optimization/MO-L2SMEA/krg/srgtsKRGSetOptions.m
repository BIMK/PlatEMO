function srgtOPT = srgtsKRGSetOptions(P, T, FIT_Fn, ...
   KRG_RegressionModel, KRG_CorrelationModel, KRG_Theta0, KRG_LowerBound, KRG_UpperBound)
%Function srgtsKRGSetOptions creates the SURROGATES Toolbox option
%structure for kriging models. This structure contains the following fiels:
%
%* GENERAL PARAMETERS
%
%   SRGT   - Identifier of the surrogate technique: 'KRG'.
%   P      - NPOINTS-by-NDV matrix, where NPOINTS is the number of points
%            of the sample and NDV is the number of design variables.
%            Default: Empty matrix.
%   T      - NPOINTS-by-1 vector of responses on the P matrix points.
%            Default: Empty vector.
%   FIT_Fn - Function handle of the fitting function (which is used to
%            optimize KRG_Theta). [@dace_fit | @srgtsXVFit].
%            Default: @dace_fit.
%
%* KRIGING PARAMETERS
%
%   KRG_RegressionModel  - Function handle to a regression model. [
%                          function_handle | @dace_regpoly0 |
%                          @dace_regpoly1 | @dace_regpoly2]. Default:
%                          @dace_regpoly0.
%   KRG_CorrelationModel - Function handle to a correlation model. [
%                          function_handle | @dace_corrcubic |
%                          @dace_correxp | @dace_correxpg |
%                          @dace_corrgauss | @dace_corrlin |
%                          @dace_corrspherical | @dace_corrspline ].
%                          Default: @dace_corrgauss.
%   KRG_Theta0           - Initial guess for theta (correlation function
%                          parameters). Default:
%                          (NPOINTS^(-1/NDV))*ones(1, NDV).
%   KRG_LowerBound       - Lower bound for theta. Default: Empty vector.
%   KRG_UpperBound       - Upper bound for theta. Default: Empty vector.
%
%The SURROGATES Toolbox uses the DACE toolbox of Lophaven et al. (2002) to
%execute the kriging algorithm. As in DACE, when KRG_LowerBound and
%KRG_UpperBound are empty ([]) there will be NO optimization on theta
%(correlation function parameters).
%
%This is how you can use srgtsKRGSetOptions:
%
%     OPTIONS = srgtsKRGSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsKRGSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsKRGSetOptions(P, T, ..
%     KRG_UpperBound, FIT_Fn, FIT_LossFn KRG_RegressionModel, ...
%     KRG_CorrelationModel, KRG_Theta0, KRG_LowerBound):
%     it creates a structure with each of the specified fields.
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
%     options = srgtsKRGSetOptions(X, Y)
%
%     options =
%
%                     SRGT: 'KRG'
%                        P: [5x1 double]
%                        T: [5x1 double]
%      KRG_RegressionModel: @dace_regpoly0
%     KRG_CorrelationModel: @dace_corrgauss
%               KRG_Theta0: 0.2000
%           KRG_LowerBound: []
%           KRG_UpperBound: []
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
srgtOPT.SRGT = 'KRG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options
switch nargin
    case 0
        srgtOPT.P = [];
        srgtOPT.T = [];
        
        srgtOPT.FIT_Fn = [];
        
        srgtOPT.KRG_RegressionModel  = [];
        srgtOPT.KRG_CorrelationModel = [];
        srgtOPT.KRG_Theta0           = [];
        srgtOPT.KRG_LowerBound       = [];
        srgtOPT.KRG_UpperBound       = [];
        
    case 2
        [npoints nvariables] = size(P);
        
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn = @dace_fit;
        
        srgtOPT.KRG_RegressionModel  = @dace_regpoly0;
        srgtOPT.KRG_CorrelationModel = @dace_corrgauss;
        srgtOPT.KRG_Theta0           = (npoints^(-1/nvariables))*ones(1, nvariables);
        srgtOPT.KRG_LowerBound       = [];
        srgtOPT.KRG_UpperBound       = [];

    otherwise
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn = FIT_Fn;
        
        srgtOPT.KRG_RegressionModel  = KRG_RegressionModel;
        srgtOPT.KRG_CorrelationModel = KRG_CorrelationModel;
        srgtOPT.KRG_Theta0           = KRG_Theta0;
        srgtOPT.KRG_LowerBound       = KRG_LowerBound;
        srgtOPT.KRG_UpperBound       = KRG_UpperBound;
end

return

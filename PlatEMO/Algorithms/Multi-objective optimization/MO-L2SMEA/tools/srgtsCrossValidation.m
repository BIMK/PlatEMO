function [PRESSRMS, eXV, yhatXV, predvarXV, srgtSRGTXV, srgtOPTXV] = srgtsCrossValidation(srgtOPT, kfolds, idxFolds)
%Function srgtsCrossValidation computes cross-validation measures for the
%given surrogate models. Thus, for example:
%
%     PRESSRMS = srgtsCrossValidation(OPTIONS): returns a 1 x 1 vector of
%     PRESSRMS value for each surrogate, which is given by:
%
%     PRESSRMS = sqrt ( PRESS / NBPOINTS );
%
%     where:
%         * PRESS (Prediction Sum of Squares): PRESS = (eXV.')*eXV
%         * OPTIONS is the SURROGATES Toolbox option structure given by
%           srgtsOptionSet, and SURROGATE is the respective surrogate
%           structure.
%
%     PRESSRMS = srgtsCrossValidation(OPTIONS, KFOLDS, IDXFOLDS):
%     does the computation of the cross-validation errors using k-fold
%     strategy. KFOLDS is the number of folds and IDXFOLDS is the vector of
%     indexes organized for NBPOINTS/k clusters. IDXFOLDS is obtained using
%     srgtsGetKfolds.
%
%     [PRESSRMS, eXV] = srgtsCrossValidation(...):
%     also returns a NBPOINTS-by-NBSURROGATES  matrix with the
%     cross-validation errors. Here the cross-validation errors obtained
%     when one data point is ignored and the surrogate is fitted to the
%     other (NBPOINTS - 1) points, with the procedure repeated for each
%     data point. In a point, the cross-validation error is:
%
%     eXV = Yhat - Y;
%
%     where Yhat is the value evaluated with the surrogate and Y is the
%     actual value of the function.
%
%     [PRESSRMS, eXV, yhatXV] = srgtsCrossValidation(...):
%     also returns the NBPOINTS-by-NBSURROGATES matrix of cross-validation
%     prediction.
%
%     [PRESSRMS, eXV, yhatXV, predvarXV] = srgtsCrossValidation(...):
%     also returns the NBPOINTS-by-NBSURROGATES matrix of cross-validation
%     prediction variance (prediction variance is implemented only for KRG
%     and PRS models; for all others it is going to be NaN).
%
%     [PRESSRMS, eXV, yhatXV, predvarXV, srgtSRGTXV] = srgtsCrossValidation(...):
%     also returns all SURROGATES created during cross-validation.
%
%     [PRESSRMS, eXV, yhatXV, predvarXV, srgtSRGTXV, ...
%     srgtOPTXV] = srgtsCrossValidation(...): also returns all
%     OPTION structures created during cross-validation.
%
%Example:
%     % basic information about the problem
%     myFN = @cos; % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     NbPoints = 5;
%     X = linspace(designspace(1), designspace(2), NbPoints)';
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     srgtOPT = srgtsPRSSetOptions(X, Y);
%
%     % calculate cross validation errors, and PRESSRMS
%     [PRESSRMS, eXV] = srgtsCrossValidation(srgtOPT);
%
%     PRESSRMS =
%
%     1.1465
%
%     eXV =
%
%     1.5700
%    -0.8006
%     0.6003
%    -0.8006
%     1.5700
%
% Results may change for different "srgtOPT" structure.

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
% check inputs
NbPoints = length(srgtOPT.T);
if nargin == 1
    kfolds   = NbPoints;
    idxFolds = [1:kfolds].';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run
Pbkp = srgtOPT.P;
Tbkp = srgtOPT.T;

NbPointsPerFold = NbPoints / kfolds; % which is also the number of clusters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create data structures
srgtSRGTXV = cell(kfolds, 1);
srgtOPTXV  = cell(kfolds, 1);
predvarXV  = NaN(NbPoints, 1);
yhatXV     = zeros(NbPoints, 1);
for c1 = 1 : kfolds
    idx = idxFolds( [ ((c1 - 1)*NbPointsPerFold + 1) : c1*NbPointsPerFold ] );

    Ptest = Pbkp(idx,:);

    Ptraining = Pbkp; Ptraining(idx,:) = [];
    Ttraining = Tbkp; Ttraining(idx)   = [];

    % new surrogates
    srgtOPT.P = Ptraining;
    srgtOPT.T = Ttraining;
    eval(sprintf('srgtSRGT = srgts%sFit(srgtOPT);', srgtOPT.SRGT));

    srgtSRGTXV{c1} = srgtSRGT;
    srgtOPTXV{c1}  = srgtOPT;

    eval(sprintf('yhatAux = srgts%sEvaluate(Ptest, srgtSRGT);', srgtOPT.SRGT));
    yhatXV(idx,:)    = yhatAux;

    switch srgtOPT.SRGT
        case {'KRG' 'GP'}
            eval(sprintf('predvarAux = srgts%sPredictionVariance(Ptest, srgtSRGT);', srgtOPT.SRGT));
            predvarXV(idx,:) = predvarAux;
        case 'PRS'
            eval(sprintf('predvarAux = srgts%sPredictionVariance(Ptest, Ptraining, srgtSRGT);', srgtOPT.SRGT));
            predvarXV(idx,:) = predvarAux;
    end
end

% cross-validation errors and PRESSRMS
eXV      = yhatXV - Tbkp;
PRESSRMS = sqrt(mean(eXV.^2));

return

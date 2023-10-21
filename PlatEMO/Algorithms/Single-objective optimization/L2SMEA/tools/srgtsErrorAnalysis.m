function CRITERIA = srgtsErrorAnalysis(srgtOPT, srgtSRGT, varargin)
%Function srgtsErrorAnalysis calculates some metrics associated with the
%quality of the fit and prediction capability of the surrogate. Thus, for
%example:
%
%     CRITERIA = srgtsErrorAnalysis(srgtOPT, srgtSRGT, Y, YHAT):
%     returns CRITERIA, a structure with the following information:
%
%          * source   : whether data or test points were used.
%          * npoints : number of points.
%          * Ystd     : Standard deviation of data.
%          * Yrange   : Range of data.
%          * RMSE     : Root mean square error.
%          * MinE     : Minimum error.
%          * MaxE     : Maximum error.
%          * MAE      : Maximum absolute error.
%          * RMSE     : Relative root mean square error.
%          * MinE     : Relative minimum error.
%          * MaxE     : Relative maximum error.
%          * MAE      : Relative maximum absolute error.
%
%          Fields meaningful only for polynomial response surfaces:
%          * PRS_SE          : Standard error.
%          * PRS_R2          : Coefficient of determination (R-square).
%          * PRS_R2a         : Adjusted R-square.
%
%          Fields meaningful only for kriging:
%          * KRG_ProcVar : Kriging process variance.
%
%     CRITERIA = srgtsErrorAnalysis(srgtOPT, srgtSRGT, Y, YHAT, X, DESIGNSPACE):
%     returns the same CRITERIA structure. However it corrects the
%     calculation of the RMSE to take into account the effect of
%     trapezoidal integration. X is the grid sample and DESIGNSPACE is the
%     matrix that defines the used to sample X.
%
%Example:
%     % basic information about the problem
%     myFN = @cos; % this could be any user-defined function
%     DesignSpace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     NbVariables = length(DesignSpace(1,:));
%
%     % create DOE
%     npoints = 5;
%     X = linspace(DesignSpace(1), DesignSpace(2), npoints)';
%     Y = feval(myFN, X);
%
%     % create test points
%     npointstest = 20;
%     Xtest = linspace(DesignSpace(1), DesignSpace(2), npointstest)';
%     Ytest = feval(myFN, Xtest);
%
%     % fit srgtSRGT model
%     srgtOPT  = srgtsPRSSetOptions(X, Y);
%     srgtSRGT = srgtsPRSFit(srgtOPT);
%
%     % check performance with data points
%     YhatData = srgtsPRSEvaluate(X, srgtSRGT);
%
%     ErrorAnalysisData = srgtsErrorAnalysis(srgtOPT, srgtSRGT, Y, YhatData)
%
%     ErrorAnalysisData =
%
%            npoints: 5
%                Ystd: 0.8367
%              Yrange: 2
%                RMSE: 0.2138
%                MinE: -0.2286
%                MaxE: 0.3429
%                 MAE: 0.3429
%               RRMSE: 0.2556
%               RMinE: -0.2732
%               RMaxE: 0.4098
%                RMAE: 0.4098
%              RYYhat: 0.9583
%              PRS_SE: 0.3381
%              PRS_R2: 0.9184
%             PRS_R2a: 0.8367
%         KRG_ProcVar: NaN
%
%     % check performance with test points
%     YhatTest = srgtsPRSEvaluate(Xtest, srgtSRGT);
%
%     ErrorAnalysisTest = srgtsErrorAnalysis(srgtOPT, srgtSRGT, Ytest, YhatTest)
%
%     ErrorAnalysisTest =
%
%            npoints: 20
%                Ystd: 0.7416
%              Yrange: 1.9864
%                RMSE: 0.2657
%                MinE: -0.4016
%                MaxE: 0.3340
%                 MAE: 0.4016
%               RRMSE: 0.3582
%               RMinE: -0.5415
%               RMaxE: 0.4503
%                RMAE: 0.5415
%              RYYhat: 0.9518
%              PRS_SE: 0.3381
%              PRS_R2: 0.8649
%             PRS_R2a: 0.8490
%         KRG_ProcVar: NaN

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
% define auxiliar variables
Y    = varargin{1};
Yhat = varargin{2};

npoints = length(Y);
w        = ones(npoints,1);
switch nargin
    case 4 % srgtOPT, srgtSRGT, Y, Yhat
        NbVariables = length(srgtOPT.P(1,:));
        X = srgtOPT.P;
    case 5 % srgtOPT, srgtSRGT, Y, Yhat, X
        X = varargin{3};
        NbVariables = length(X(1,:));
    case 6 % srgtOPT, srgtSRGT, Y, Yhat, X, DesignSpace
        X = varargin{3};
        NbVariables = length(X(1,:));
        Xaux = srgtsScaleVariable(X, ...
            varargin{4}, ...
            [zeros(1,NbVariables); ones(1,NbVariables)]);
        
        npoints = length(X(:,1));
        w = ones(size(npoints,1));
        for c1 = 1 : npoints % correction for trapezoidal integration
            x = Xaux(c1,:);
            for c2 = 1 : NbVariables
                if ( ( x(c2) == 0 ) || ( x(c2) == 1 ) )
                    w(c1) = 0.5*w(c1);
                end
            end
        end
end

errors = Yhat - Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign outputs
%        npoints: 'Number of points'
%            Ystd: 'Standard deviation of data'
%          Yrange: 'Range of data'
%            RMSE: 'Root mean square error'
%            MinE: 'Minimum error'
%            MaxE: 'Maximum error'
%             MAE: 'Maximum absolute error'
%           RRMSE: 'Relative root mean square error'
%           RMinE: 'Relative minimum error'
%           RMaxE: 'Relative maximum error'
%            RMAE: 'Relative maximum absolute error'
%          RYYhat: 'Correlation coefficient between Y and Yhat'


CRITERIA.npoints = npoints;

CRITERIA.Ystd   = std(Y);     % Standard deviation of data
CRITERIA.Yrange = range(Y); % Range of data

CRITERIA.RMSE = sqrt( mean( w.*(errors.^2) ) );
CRITERIA.MinE = min(errors);
CRITERIA.MaxE = max(errors);
CRITERIA.MAE  = max( abs(errors) );

CRITERIA.RRMSE = CRITERIA.RMSE/CRITERIA.Ystd;
CRITERIA.RMinE = CRITERIA.MinE/CRITERIA.Ystd;
CRITERIA.RMaxE = CRITERIA.MaxE/CRITERIA.Ystd;
CRITERIA.RMAE  = CRITERIA.MAE/CRITERIA.Ystd;

% Criteria for individual srgtSRGTs (basically, polynomial response
% surfaces and kriging)
% For polynomial response surfaces:
%              SE: 'Standard error'
%              R2: 'R2 (coefficient of determination)'
%             R2a: 'R2a (adjusted R2)'
%
% For kriging:
%           KRGPV: 'Kriging process variance'

CRITERIA.PRS_SE  = NaN;
CRITERIA.PRS_R2  = NaN;
CRITERIA.PRS_R2a = NaN;
CRITERIA.KRG_ProcVar = NaN;

switch srgtOPT.SRGT
    case 'PRS'
        MeanY   = mean(Y); % Mean of the response.
        SST     = sum( (Y - MeanY).^2 ); % Total sum of squares
        NbCoeff = length(srgtSRGT.PRS_Beta);
        
        CRITERIA.PRS_SE = srgtSRGT.PRS_SE;
        
        X   = srgtsPRSCreateGramianMatrix(X, NbVariables, srgtOPT.PRS_Degree, srgtSRGT.PRS_RemovedIdx);
        SSE = sum( errors.^2 );          % Sum of squared errors
        SSR = sum(abs(Yhat - MeanY).^2); % Regression sum of squares
        
        CRITERIA.PRS_R2  = (SST - SSE)/SST; % Coefficient of determination (R-square)
        CRITERIA.PRS_R2a = 1 - ( (npoints - 1) ./ (npoints - NbCoeff) ).*(1 - CRITERIA.PRS_R2); % adjusted R-square
        
    case 'KRG'
        CRITERIA.KRG_ProcVar = srgtSRGT.KRG_DACEModel.sigma2;
end

return

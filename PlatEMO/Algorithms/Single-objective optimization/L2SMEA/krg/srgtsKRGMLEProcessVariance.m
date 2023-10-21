function MLEPV = srgtsKRGMLEProcessVariance(srgtOPT, srgtSRGT)
%Function srgtsKRGMLEProcessVariance returns the maximum likelihood
%estimate of the process variance.
%
%    MLEPV = srgtsKRGMLEProcessVariance(srgtOPT, srgtSRGT)
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
%     % option for kriging
%     srgtOPT  = srgtsKRGSetOptions(X, Y);
%     srgtSRGT = srgtsKRGFit(srgtOPT);
% 
%     MLEPV = srgtsKRGMLEProcessVariance(srgtOPT, srgtSRGT)
% 
%     MLEPV =
% 
%     11.4553

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

[dummy01, dummy02, MLEPV] = dace_fit(...
        srgtOPT.P, ...
        srgtOPT.T, ...
        srgtOPT.KRG_RegressionModel, ...
        srgtOPT.KRG_CorrelationModel, ...
        srgtSRGT.KRG_DACEModel.theta);


return

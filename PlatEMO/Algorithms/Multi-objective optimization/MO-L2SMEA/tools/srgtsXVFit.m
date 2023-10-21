function [srgtSRGT srgtSTT] = srgtsXVFit(srgtOPT)
%Function srgtsXVFit can be used to fit Gaussian process, kriging, and
%radial basis function models using cross validation. Please, check the
%options' structure of these surrogates for details.
%
%
%Although this function was intended to be used by the fitting functions of
%Gaussian process, kriging, and radial basis function models, you can also
%use it as
%
%    [srgtSRGT srgtSTT] = srgtsXVFit(srgtOPT)
%
%srgtSRGT is the surrogate structure and srgtSTT is a structure that
%contains the following fields:
%* FIT_FnVal: PRESSRMS
%* For kriging models:
%    * KRG_DACEPerf (see srgtsKRGFit for details)


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

srgtSTT = srgtsFitCreateState(srgtOPT);

eval(sprintf('lowerbound = srgtOPT.%s_LowerBound;', srgtOPT.SRGT));
eval(sprintf('upperbound = srgtOPT.%s_UpperBound;', srgtOPT.SRGT));
eval(sprintf('srgtOPT.%s_LowerBound = [];', srgtOPT.SRGT));
eval(sprintf('srgtOPT.%s_UpperBound = [];', srgtOPT.SRGT));

switch srgtOPT.SRGT
    case 'KRG'
        srgtOPT.FIT_Fn = @dace_fit;
    case 'RBF'
        srgtOPT.FIT_Fn = @rbf_build;
    case 'GP'
        srgtOPT.FIT_Fn = @gpml_gpr_fit;
end

ndv = length(lowerbound);
if ndv < 10
    popsize = max(20*ndv, 50);
else
    popsize = 100;
end

foptm   = Inf;
itermax = 100; % maximum of function evaluations = itermax*popsize
for c1 = 1 : 4
    [ptemp ftemp] = srgtsOPTMDE(@srgtsXVFittingObjective, -Inf, ndv, ...
        lowerbound, upperbound, srgtOPT, popsize, itermax, ...
        0.8, 0.8, 7, 0); % DE parameters
    if ftemp < foptm
        foptm = ftemp;
        poptm = ptemp;
    end
end

switch srgtOPT.SRGT % create the srgtSRGT structure
    case 'GP'
        srgtOPT.GP_LogTheta0 = poptm';
        [srgtSRGT srgtSTT] = srgtsGPFit(srgtOPT);

    case 'KRG'
        srgtOPT.KRG_Theta0 = poptm;
        srgtSRGT = srgtsKRGFit(srgtOPT);
        srgtSTT.KRG_DACEPerf.nv   = [];
        srgtSTT.KRG_DACEPerf.perf = [];
        
    case 'RBF'
        srgtOPT.RBF_c = poptm;
        srgtSRGT = srgtsRBFFit(srgtOPT);

end

srgtSTT.FIT_FnVal = sqrt(foptm);

return

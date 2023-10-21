function PRESSRMS = srgtsXVFittingObjective(x, srgtOPT)
% not meaningful for the user

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

switch srgtOPT.SRGT % prepare for evaluation of the loss function
    case 'KRG'
        srgtOPT.KRG_Theta0 = x;
    case 'RBF'
        srgtOPT.RBF_c = x;
    case 'GP'
        srgtOPT.GP_LogTheta0 = x';
end

PRESSRMS = srgtsCrossValidation(srgtOPT);

PRESS    = PRESSRMS.^2;

return

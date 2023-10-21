function kfolds = srgtsGetKfoldValues(NbPoints)
%Function srgtsGetKfoldValues calculates the possible kfold values for a
%certain number of points. Thus, for example:
%
%     KFOLDS = srgtsGetKfoldValues(NPOINTS): returns a two column matrix
%     in which the first column is the kfold value and the second column is
%     the number of points in each fold.
%     It is important to notice that k = 1 is not an option.
%
%Example:
%     NPOINTS = 200;
%     KFOLDS  = srgtsGetKfoldValues(NPOINTS)
%
%     KFOLDS =
% 
%           2   100
%           4    50
%           5    40
%           8    25
%           10   20
%           20   10
%           25    8
%           40    5
%           50    4
%           100   2
%           200   1
%
%See also srgtsGetKfolds and srgtsComputeCrossValidation.

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
% run
kfolds = [];
for c1 = 2 : NbPoints
    candidate = NbPoints/c1;
    if round(candidate) == candidate
        kfolds = vertcat(kfolds, [c1 candidate]);
    end
    
end

return

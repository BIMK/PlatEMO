function IdxFolds = srgtsGetKfoldMaxMin(P, NbPoints, kfold, NbPointsFold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IdxFolds = [];

NbDone    = 0;
RemDoE    = P;
NbRemPoints = NbPoints;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculations
while NbPoints > NbDone,

    [xfold, dummy01, RemDoE, RPIdx] = srgtsDOESubSample(RemDoE, ...
        NbPointsFold, ...
        'MaxMin', 1000);

    NbRemPoints = length(RPIdx);

    for c1 = 1 : NbPointsFold
        for c2 = 1 : NbPoints
            if isequal( xfold(c1,:), P(c2,:) )
                IdxFolds = vertcat(IdxFolds, c2);
            end
        end
    end

    NbDone = NbDone + NbPointsFold;

end

return

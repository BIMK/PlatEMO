function Population = EnvironmentalSelection(Population,R,N)
% The environmental selection of RPEA

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the Tchebychev distances
    Distance = TchebychevDistance(Population.objs,R);
    
    %% Environmental selection
    RemainP = 1 : length(Population);
    RemainR = 1 : size(R,1);
    while length(RemainP) > length(Population)-N
        if isempty(RemainR)
            RemainR = 1 : size(R,1);
        end
        [temp,imin] = min(Distance(RemainP,RemainR),[],1);
        [~,jmin]    = min(temp);
        imin        = imin(jmin);
        RemainP(imin) = [];
        RemainR(jmin) = [];
    end
    Population = Population(setdiff(1:length(Population),RemainP));
end
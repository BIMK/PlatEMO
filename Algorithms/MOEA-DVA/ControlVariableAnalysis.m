function [DiverIndexes,ConverIndexes] = ControlVariableAnalysis(Global,NCA)
% Control variable analysis

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    DiverIndexes  = false(1,Global.D);
    ConverIndexes = false(1,Global.D);
    for i = 1 : Global.D
        x      = unifrnd(Global.lower,Global.upper);
        S      = repmat(x,NCA,1);
        S(:,i) = ((1:NCA)'-1+rand(NCA,1))/NCA*(Global.upper(i)-Global.lower(i)) + Global.lower(i);
        S      = INDIVIDUAL(S);
        [~,MaxFNo] = NDSort(S.objs,inf);
        if MaxFNo == length(S)
            ConverIndexes(i) = true;
        else
            DiverIndexes(i) = true;
        end
    end
end
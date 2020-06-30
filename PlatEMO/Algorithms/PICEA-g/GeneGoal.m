function NewGoal = GeneGoal(PopObj,NGoal)
% Generate new goals

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Gmax    = max(PopObj,[],1)*1.2;
    Gmin    = min(PopObj,[],1);
    NewGoal = unifrnd(repmat(Gmin,NGoal,1),repmat(Gmax,NGoal,1));
end
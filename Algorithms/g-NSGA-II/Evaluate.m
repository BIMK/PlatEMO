function PopObj = Evaluate(PopObj,Point)
% g-dominance based function evaluation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Point = repmat(Point,size(PopObj,1),1);
    Flag  = all(PopObj<=Point,2) | all(PopObj>=Point,2);
    Flag  = repmat(Flag,1,size(PopObj,2));
    PopObj(~Flag) = PopObj(~Flag) + 1e10;
end
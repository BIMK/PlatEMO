function [Population,Dec,Mask] = Prediction(Problem,ChangeCount,DecSource,MaskSource)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    P = 4;
    if ChangeCount < P
        Mask = MaskSource{ChangeCount+1};
    else
        Mask = SVR(MaskSource,ChangeCount,P);
    end
    score      = sum(Mask,1);
    Dec        = MLP(DecSource,ChangeCount,P,score);
    Population = Problem.Evaluation(Dec.*Mask);
end
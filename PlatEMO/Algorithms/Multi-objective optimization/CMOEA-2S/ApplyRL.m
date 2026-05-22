function [Prob1,Prob2] = ApplyRL(State, Model, Paras)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    TestX1 = [State,1]; % 1-th action
    TestX2 = [State,2]; % 2-th action
    ps     = Paras.ps;
    qs     = Paras.qs;
    TestX1 = mapminmax('apply',TestX1',ps);
    TestX1 = TestX1';
    TestX2 = mapminmax('apply',TestX2',ps);
    TestX2 = TestX2';

    Prob1 = testNet(TestX1,Model,Paras);
    Prob1 = mapminmax('reverse',Prob1',qs);
    Prob1 = Prob1';

    Prob2 = testNet(TestX2,Model,Paras);
    Prob2 = mapminmax('reverse',Prob2',qs);
    Prob2 = Prob2';
end
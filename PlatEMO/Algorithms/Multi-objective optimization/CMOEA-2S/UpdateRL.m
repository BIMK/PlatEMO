function [Model, Paras] = UpdateRL(Records, Model, Paras, Gama, StateNum)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % update model here
    qs     = Paras.qs;
    TrainX = Records(:,1:StateNum+1);
    % Normalization
    [TrainX,ps] = mapminmax(TrainX');
    TrainX      = TrainX';

    Prob = testNet(TrainX,Model,Paras);
    Prob = mapminmax('reverse',Prob',qs);
    Prob = Prob';
    Prob = Prob(:,1);

    % update model here
    TrainY      = Records(:,StateNum+2)+Gama*max(Prob);
    [TrainY,qs] = mapminmax(TrainY');
    TrainY      = TrainY';
    Paras.ps    = ps;
    Paras.qs    = qs;
    Model       = updatemodel(TrainX,TrainY,Paras,Model);
end
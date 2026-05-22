function [model,Paras] = BuildRL(Records,StateNum)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Build DQN model
    TrainX = Records(:,1:StateNum+1);
    TrainY = Records(:,StateNum+2:end);
    % Normalization
    [TrainX,ps]   = mapminmax(TrainX');
    TrainX        = TrainX';
    [TrainY,qs]   = mapminmax(TrainY');
    TrainY        = TrainY';
    Paras.ps      = ps;
    Paras.qs      = qs;
    [model,Paras] = trainmodel(TrainX,TrainY,Paras);
end
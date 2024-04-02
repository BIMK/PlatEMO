function [netset,psset,qsset] = TrainModel(cycle,NumEsp,trainDecs,trainLabel,Index)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yingwei Li

    subTrain = cell(1,NumEsp);
    netset   = cell(1,NumEsp);
    psset    = cell(1,NumEsp);
    qsset    = cell(1,NumEsp);
    for j = 1 : cycle
        for i = 1 : NumEsp
            subTrain{i} = trainDecs(:,Index{j}==i);
            Input       = subTrain{i};
            [Input,ps]  = mapminmax(Input');
            Output      = trainLabel; 
            [Output,qs] = mapminmax(Output');
            net = feedforwardnet(10);
            net.trainParam.showWindow = false;
            net = train(net, Input, Output);
            netset{i} = net;
            psset{i}  = ps;
            qsset{i}  = qs;
        end
    end
end
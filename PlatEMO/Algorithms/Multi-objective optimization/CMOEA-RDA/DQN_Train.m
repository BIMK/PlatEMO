function net = DQN_Train(XTrain,YTrain)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    layers = [featureInputLayer(101, 'Name', 'input')
              fullyConnectedLayer(30, 'Name', 'fc1')
              reluLayer('Name', 'relu1')
              fullyConnectedLayer(1, 'Name', 'fc3')
              tanhLayer('Name', 'tanh1')];  % Define network structure
    
    net           = dlnetwork(layers);
    numEpochs     = 100;
    miniBatchSize = 50;
    learnRate     = 0.001;
    averageGrad   = [];
    averageSqGrad = [];
    
    XTrain = XTrain';
    YTrain = YTrain';
   
    for epoch = 1 : numEpochs       % train network
        for i = 1 : miniBatchSize:size(XTrain, 2)
            idx = i : min(i+miniBatchSize-1, size(XTrain, 2));
            X   = XTrain(:,idx);
            Y   = YTrain(:,idx);

            dlX     = dlarray(X, 'CB');
            dlYPred = forward(net, dlX);

           
            dlY  = dlarray(Y, 'CB');
            loss = mse(dlYPred, dlY);
    
            [~, gradients] = dlfeval(@modelGradients, net, loss);
    
            [net, averageGrad, averageSqGrad] = adamupdate(net, gradients, ...
                averageGrad, averageSqGrad, epoch, learnRate);
        end
    
    end
end

function [loss, gradients] = modelGradients(net, loss)
    gradients = dlgradient(loss, net.Learnables);
end
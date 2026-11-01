function net = DQN_Update(net,XTrain,YTrain)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    numEpochs     = 100;
    miniBatchSize = 50;
    learnRate     = 0.001;
    
    averageGrad   = [];
    averageSqGrad = [];
    
    XTrain = XTrain';
    YTrain = YTrain';

    for epoch = 1 : numEpochs   % update network
        for i = 1 : miniBatchSize:size(XTrain, 2)
            idx = i : min(i+miniBatchSize-1, size(XTrain, 2));
            X   = XTrain(:,idx);
            Y   = YTrain(:,idx);
            

            dlX = dlarray(X, 'CB');
            dlY = dlarray(Y, 'CB');

            [~, gradients] = dlfeval(@modelGradients, net, dlX, dlY);

            [net, averageGrad, averageSqGrad] = adamupdate(net, gradients, ...
                averageGrad, averageSqGrad, epoch, learnRate);
        end
    end
end


function [loss, gradients] = modelGradients(net, dlX, dlY)
    dlYPred   = forward(net, dlX);
    loss      = customLoss(dlYPred, dlY);
    gradients = dlgradient(loss, net.Learnables);
end

function loss = customLoss(Y, T)
    loss = sum((Y - T).^2, 'all');
end
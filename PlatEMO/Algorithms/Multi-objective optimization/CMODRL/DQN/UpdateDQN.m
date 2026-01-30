function net = UpdateDQN(net,XTrain,YTrain)

    numEpochs = 100;
    miniBatchSize = 60;
    learnRate = 0.001;
    
    averageGrad = [];
    averageSqGrad = [];
    
    XTrain = XTrain';
    YTrain = YTrain';
    
    for epoch = 1:numEpochs
        for i = 1:miniBatchSize:size(XTrain, 2)
            idx = i:min(i+miniBatchSize-1, size(XTrain, 2));
            X = XTrain(:,idx);
            Y = YTrain(:,idx);
            dlX = dlarray(X, 'CB');
            dlY = dlarray(Y, 'CB');
    
            [~, gradients] = dlfeval(@modelGradients, net, dlX, dlY);
    
            [net, averageGrad, averageSqGrad] = adamupdate(net, gradients, ...
                averageGrad, averageSqGrad, epoch, learnRate);
        end
    
    end
end

function [loss, gradients] = modelGradients(net, dlX, dlY)
    dlYPred = forward(net, dlX);
    loss = customLoss(dlYPred, dlY);
    gradients = dlgradient(loss, net.Learnables);
end

function loss = customLoss(Y, T)
    loss = sum((Y - T).^2, 'all');
end
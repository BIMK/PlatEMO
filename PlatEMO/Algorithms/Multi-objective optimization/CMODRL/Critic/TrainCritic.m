function net = TrainCritic(XTrain,YTrain)

    layers = [
        featureInputLayer(101, 'Name', 'input')
        fullyConnectedLayer(30, 'Name', 'fc1')
        reluLayer('Name', 'relu1')
        fullyConnectedLayer(30, 'Name', 'fc2')
        reluLayer('Name', 'relu2')
        fullyConnectedLayer(1, 'Name', 'fc3')
        tanhLayer('Name', 'tanh1')
        ];
    
    net = dlnetwork(layers);
    
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
            dlYPred = forward(net, dlX);
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

function loss = customLoss(Y, T)
    loss = sum((Y - T).^2, 'all');
end
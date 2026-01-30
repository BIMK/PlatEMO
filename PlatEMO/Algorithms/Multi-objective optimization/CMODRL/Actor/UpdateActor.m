function net = UpdateActor(net,critic_net,XTrain)

    numEpochs = 100;
    miniBatchSize = 60;
    learnRate = 0.001;
    
    averageGrad = [];
    averageSqGrad = [];
    
    XTrain = XTrain';
    
    for epoch = 1:numEpochs
        for i = 1:miniBatchSize:size(XTrain, 2)
            idx = i:min(i+miniBatchSize-1, size(XTrain, 2));
            X = XTrain(:,idx);
            dlX = dlarray(X, 'CB');
            action = forward(net,dlX);
            reward = forward(critic_net,dlarray([X;double(extractdata(action))],'CB'));
            loss = sum(-reward(:,1));
    
            [~, gradients] = dlfeval(@modelGradients1, net, loss);
    
            [net, averageGrad, averageSqGrad] = adamupdate(net, gradients, ...
                averageGrad, averageSqGrad, epoch, learnRate);
        end
    
    end
end

function [loss, gradients] = modelGradients(net, X, dlX, critic_net)
    action = forward(net, dlX);
    dlYPred = critic_net([X,double(extractdata(action))]);
    loss = customLoss(dlYPred(:,end));
    gradients = dlgradient(single(loss), net.Learnables);
end

function [loss, gradients] = modelGradients1(net, loss)
    gradients = dlgradient(loss, net.Learnables);
end

function loss = customLoss(Y)
    loss = sum(dlarray(-Y,'CB'));
end
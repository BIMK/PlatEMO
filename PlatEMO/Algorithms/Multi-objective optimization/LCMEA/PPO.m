classdef PPO < handle
% Proximal policy optimization

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        % Actor and Critic model
        CriticNetwork;
        ActorNetwork;

        % Environmental parameters
        NumAction;
        NumState;

        % Memory pool
        StateBuffer;
        ActionBuffer;
        LogProbBuffer;
        RewardBuffer;
        DoneBuffer;
        ValueBuffer;

        % Step counter used for updating model
        StepStamp = 1;
        Flag      = false;
        MaxGen    = 1;
        Gen       = 1;

        % Loss buffer
        acLoss = [];
        crLoss = [];
        reward = [];
    end

    properties(Access=private)
        % Model parameter
        ActorLearnRate  = 1e-2;
        CriticLearnRate = 2e-2;
        NumEnvs         = 1;
        NumSteps        = 20;
        BatchSize       = 20;
        MiniBatchSize   = 4;
        Gamma           = 0.99;
        GAELambda       = 0.95;
        UpdateEpochs    = 2;
        ClipCoef        = 0.1;
        EntCoef         = 0.01;
        VFCoef          = 0.5;
        MaxGradNorm     = 0.5;
        % Optimizer parameter (the SGDM solver)
        Momentum = 0.9;
    end
    methods
        function obj = PPO(NumAction, NumState, MaxGen)
            obj.MaxGen = MaxGen;

            %% Environmental parameter
            obj.NumAction = NumAction;
            obj.NumState  = NumState;

            %% Define Critic model
            criticPath = [featureInputLayer(obj.NumState,'Name','CriticInput')
                batchNormalizationLayer
                tanhLayer
                fullyConnectedLayer(60,'Name','CriticHidden', 'WeightsInitializer',@(sz) randn(sz) * sqrt(2))
                batchNormalizationLayer
                tanhLayer
                fullyConnectedLayer(1,'Name','CriticOutput', 'WeightsInitializer',@(sz) randn(sz) * sqrt(1))
                batchNormalizationLayer
                tanhLayer];            
            obj.CriticNetwork = dlnetwork(criticPath,Initialize=true);

            %% Define Actor modhel
            actorPath = [featureInputLayer(obj.NumState,'Name','ActorInput')
                batchNormalizationLayer
                tanhLayer
                fullyConnectedLayer(60,'Name','ActorHidden', 'WeightsInitializer',@(sz) randn(sz) * sqrt(2))
                batchNormalizationLayer
                tanhLayer
                fullyConnectedLayer(sum(obj.NumAction),'Name','ActorOutput', 'WeightsInitializer',@(sz) randn(sz) * sqrt(0.01))
                batchNormalizationLayer
                softmaxLayer];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
            obj.ActorNetwork = dlnetwork(actorPath,Initialize=true);

            %% Memory pool initialization
            obj.StateBuffer   = dlarray(zeros(obj.NumState, obj.NumSteps),"CB");
            obj.ActionBuffer  = dlarray(zeros(sum(obj.NumAction), obj.NumSteps),"CB");
            obj.LogProbBuffer = dlarray(zeros(1, obj.NumSteps),"CB");
            obj.RewardBuffer  = dlarray(zeros(1, obj.NumSteps),"CB");
            obj.ValueBuffer   = dlarray(zeros(1, obj.NumSteps),"CB");
            obj.DoneBuffer    = dlarray(zeros(1, obj.NumSteps),"CB");
        end

        function [obj, action] = GetAction(obj, State, Reward)
            State = dlarray(State',"CB");

            if mod(obj.StepStamp, obj.NumSteps) == 1 && obj.Gen > obj.NumSteps
                % Reset stepstamp
                stateNext = State;
                doneNext  = 0;

                %% Update reward buffer
                obj.RewardBuffer(obj.StepStamp-1) = Reward;
                
                %% Update actor and critic
                obj = obj.UpdateModel(stateNext, doneNext);

                %% Reset stepstamp
                obj.StepStamp = 1;
            else
                %% Record reward
                if obj.StepStamp > 2
                    obj.RewardBuffer(obj.StepStamp-1) = Reward;
                end                
            end
            %% record path
            obj.StateBuffer(:, obj.StepStamp)  = State;
            obj.DoneBuffer(obj.StepStamp)     = 0;
            [action, logProb, value] = obj.GetActionValue(State);
            obj.ActionBuffer(:, obj.StepStamp)  = action;
            obj.LogProbBuffer(:, obj.StepStamp) = logProb;
            obj.ValueBuffer(obj.StepStamp)     = value;

            %% Update generation recorder
            obj.StepStamp = obj.StepStamp + 1;
            obj.Gen = obj.Gen + 1;
            obj.reward = [obj.reward,Reward];
        end

        function [action, logProb, value] = GetActionValue(obj, State)
            % Predict value
            value = obj.GetValue(State);

            % Predict action and other parameters
            logits = forward(obj.ActorNetwork, State);
            if obj.Gen == 1 || rand() < 0.2
                action = randperm(obj.NumAction, 1);
            else
                action = obj.SampleAction(logits);
            end
            logProb = obj.CalLogProbs(logits);
            logProb = sum(logProb, 1);
        end

        function value = GetValue(obj, State)
            State = dlarray(State, "CB");
            value = forward(obj.CriticNetwork, State);
        end

        function obj = UpdateModel(obj, StateNext, DoneNext)
            %% Boostrap value
            valueNext = obj.GetValue(StateNext);
            returns = zeros(size(obj.RewardBuffer));
            for i = obj.NumSteps : -1 : 1
                if i == obj.NumSteps
                    nextNonTerminal = 1 - DoneNext;
                    nextReturn = valueNext;
                else
                    nextNonTerminal = 1 - obj.DoneBuffer(i+1);
                    nextReturn = returns(i+1);
                end
                returns(i) = obj.RewardBuffer(i) + obj.Gamma*nextNonTerminal*nextReturn;
            end
            advantages = returns - obj.ValueBuffer;

            %% Optimize the actor and critic network
            frac = 1 - (obj.Gen-1) / (2*obj.MaxGen) ;
            acLearnRate = frac * obj.ActorLearnRate;
            crLearnRate = frac * obj.CriticLearnRate;
            for epoch = 1 : 1 : obj.UpdateEpochs
                randIndex  = randperm(obj.BatchSize);
                velActor   = [];
                velCritic  = [];
                for start = 1 : obj.MiniBatchSize : obj.BatchSize
                    trainIndex = randIndex(start:1:start-1+obj.MiniBatchSize);
                    trainAdvantages = advantages(trainIndex);
                    trainAdvantages = (trainAdvantages-mean(trainAdvantages)) ./ (std(trainAdvantages)+1e-8);

                    %% Calculate gradient and loss
                    [acloss, acGradients] = dlfeval(@obj.ActorLoss,obj.ActorNetwork, obj.StateBuffer(:, trainIndex), trainAdvantages, obj.LogProbBuffer(:,trainIndex));
                    %% Optimize actor using the SGDM optimizer
                    [obj.ActorNetwork, velActor] = sgdmupdate(obj.ActorNetwork, acGradients, velActor, acLearnRate, obj.Momentum);

                    %% Calculate gradient and loss
                    % crGradients = dlgradient(vLoss*obj.VFCoef,obj.CriticNetwork.Learnables);
                    [crloss, crGradients] = dlfeval(@obj.CriticLoss,obj.CriticNetwork, obj.StateBuffer(:, trainIndex), obj.ValueBuffer(trainIndex), returns(trainIndex));
                    %% Optimize critic using the SGDM optimizer
                    [obj.CriticNetwork, velCritic] = sgdmupdate(obj.CriticNetwork, crGradients, velCritic, crLearnRate, obj.Momentum);
                end
            end
        end

        function [crLoss, crGradients] = CriticLoss(obj, CriticNet, State, Values, Returns)
            value = forward(CriticNet, State);

            vLossUnclipped = (value-Returns).^2;
            vClipped       = Values + min(max(value-Values, -obj.ClipCoef), obj.ClipCoef);
            vLossClipped   = (vClipped-Returns).^2;
            vLossMax       = max(vLossClipped,vLossUnclipped);
            crLoss         = 0.5*mean(vLossMax);
            crLoss = crLoss*obj.VFCoef;

            % Calculate gradient
            crGradients = dlgradient(crLoss,CriticNet.Learnables);
        end

        function [acLoss, acGradients] = ActorLoss(obj, ActorNet, State, Advantages, LogProb)
            logits = forward(ActorNet, State);

            newLogProb = obj.CalLogProbs(logits);
            newLogProb = sum(newLogProb,1);
            entropy    = obj.CalEntropy(logits);
            entropy    = sum(entropy,1);
            logRatio = newLogProb - LogProb;
            ratio    = exp(logRatio);

            acLoss1 = -Advantages .* ratio;
            acLoss2 = -Advantages .* min(max(ratio, 1-obj.ClipCoef),1+obj.ClipCoef);
            acLoss  = mean(max(acLoss1, acLoss2)) - obj.EntCoef*mean(entropy);
            acGradients = dlgradient(acLoss, ActorNet.Learnables);
        end

        function logProbs = CalLogProbs(obj, logits)
            logProbs = logits - log(sum(exp(logits), 1));
        end

        function entropy = CalEntropy(obj, logits)
            probs = exp(logits) ./ sum(exp(logits), 1);
            entropy  = - sum(probs.*log(probs),1);
        end

       function action = SampleAction(obj, logits)
           probs = extractdata(exp(logits) ./ sum(exp(logits), 1))';
           probs(end) = 1-sum(probs(1:end-1));
           % sample
           action  = randsrc(1,1,[1:obj.NumAction;probs]);
        end
    end
end
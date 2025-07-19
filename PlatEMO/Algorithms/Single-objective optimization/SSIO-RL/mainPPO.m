function mainPPO(TrainN, TrainFE, TrainOffspring, Problem, fileName)
% Environment

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    ObservationInfo             = rlNumericSpec([6 1]);
    ObservationInfo.Name        = 'myobs';
    ObservationInfo.Description = 'Pdec';
    
    ActionInfo             = rlNumericSpec([31 1]);
    ActionInfo.LowerLimit  = repmat(-1,1,31)';
    ActionInfo.UpperLimit  = ones(1,31)';
    ActionInfo.Name        = 'Evolution Action';
    ActionInfo.Description = 'Weight';

    func_handle_myResetFunction = @() myResetFunction(TrainN, TrainFE, TrainOffspring, Problem);
    % Create environment
    env = rlFunctionEnv(ObservationInfo,ActionInfo,'myStepFunction', func_handle_myResetFunction);
    
    obsInfo         = getObservationInfo(env);
    obsInfo.Name    = 'myobs';
    numObservations = obsInfo.Dimension(1);

    actInfo         = getActionInfo(env);
    numActions      = actInfo.Dimension(1);
    
    criticNetwork = [imageInputLayer([numObservations 1],'Normalization','none','Name','myobs')
                     fullyConnectedLayer(32,'Name','CriticFC1')
                     reluLayer('Name','CriticRelu1')
                     fullyConnectedLayer(64,'Name','CriticFC2')
                     reluLayer('Name','CriticRelu2')
                     fullyConnectedLayer(128,'Name','CriticFC3')
                     reluLayer('Name','CriticRelu3')
                     fullyConnectedLayer(64,'Name','CriticFC4')
                     reluLayer('Name','CriticRelu4')
                     fullyConnectedLayer(32,'Name','CriticFC5')
                     reluLayer('Name','CriticRelu5')
                     fullyConnectedLayer(1,'Name','CriticOutput')];
    
    critic = rlValueFunction(criticNetwork, obsInfo, 'UseDevice', 'cpu');
    
    actorNetwork = layerGraph([featureInputLayer(numObservations, 'Normalization', 'none', 'Name', 'myobs')
                               fullyConnectedLayer(32, 'Name', 'fc1')
                               reluLayer('Name', 'relu1')
                               fullyConnectedLayer(64, 'Name', 'fc2')
                               reluLayer('Name', 'relu2')
                               fullyConnectedLayer(128, 'Name', 'fc3')
                               reluLayer('Name', 'relu3')
                               fullyConnectedLayer(64, 'Name', 'fc4')
                               reluLayer('Name', 'relu4')
                               fullyConnectedLayer(32, 'Name', 'fc5')
                               reluLayer('Name', 'relu5')]);
    
    meanLayer    = fullyConnectedLayer(numActions, 'Name', 'mean');
    actorNetwork = addLayers(actorNetwork, meanLayer);
    actorNetwork = connectLayers(actorNetwork, 'relu5', 'mean');
    
    tanhLayer1   = tanhLayer('Name', 'tanh');
    actorNetwork = addLayers(actorNetwork, tanhLayer1);
    actorNetwork = connectLayers(actorNetwork, 'mean', 'tanh');
    
    stdLayer     = fullyConnectedLayer(numActions, 'Name', 'std');
    actorNetwork = addLayers(actorNetwork, stdLayer);
    actorNetwork = connectLayers(actorNetwork, 'relu5', 'std');

    softplusLayer1 = softplusLayer('Name', 'softplus');
    actorNetwork   = addLayers(actorNetwork, softplusLayer1);
    actorNetwork   = connectLayers(actorNetwork, 'std', 'softplus');

    actor = rlContinuousGaussianActor(actorNetwork, obsInfo, actInfo, ...
                                      'ActionMeanOutputNames', 'tanh', ...
                                      'ActionStandardDeviationOutputNames', 'softplus', ...
                                      'ObservationInputNames', 'myobs');
        
    opt = rlPPOAgentOptions('SampleTime', 1, ...
                            'DiscountFactor', 0.99, ...
                            'EntropyLossWeight', 0.01, ...
                            'ExperienceHorizon', 512, ...
                            'MiniBatchSize', 64, ...
                            'NumEpoch', 3, ...
                            'MaxMiniBatchPerEpoch', 100, ...
                            'LearningFrequency', 64, ...
                            'ClipFactor', 0.2, ...
                            'AdvantageEstimateMethod', "gae", ...
                            'GAEFactor', 0.95, ...
                            'NormalizedAdvantageMethod', "none", ...
                            'AdvantageNormalizingWindow', 1e6);
    
    agent      = rlPPOAgent(actor,critic,opt);
    scriptPath = fileparts(mfilename('fullpath'));
    trainOpts  = rlTrainingOptions('MaxEpisodes',1000,...
                                   'MaxStepsPerEpisode',200,...
                                   'ScoreAveragingWindowLength',30,...
                                   'Verbose',true,...
                                   'Plots','training-progress',...
                                   'StopTrainingCriteria','AverageReward',...
                                   'StopTrainingValue',300,...
                                   'Plots','none');
    train(agent,env,trainOpts);
    save(fullfile(scriptPath, fileName), 'agent');
end
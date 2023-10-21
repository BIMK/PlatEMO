classdef DQN < handle
% DQN network with gamma being zero

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        % Constants
        MEMORY_CAPACITY = 512
        BATCH_SIZE = 8
        LR = 0.01             	% learning rate
        EPSILON = 0.9         	% greedy policy
        TARGET_REPLACE_ITER = 7	% target update frequency
        % Variables
        net
        n_states              	% state number
        n_actions            	% action number
        learn_step_counter      % for target updating
        memory_counter          % for storing memory
        memory_scurr            % current state
        memory_a              	% action
        memory_r             	% reward
        memory_snext          	% next state
    end
    methods
        function obj = DQN(inDim, outDim)
        % Network definition
            layers = [featureInputLayer(inDim, 'Name', 'input')
                      fullyConnectedLayer(128, 'Name', 'fc_1')
                      reluLayer('Name', 'relu_1')
                      fullyConnectedLayer(256, 'Name', 'fc_2')
                      reluLayer('Name', 'relu_2')
                      fullyConnectedLayer(128, 'Name', 'fc_3')
                      reluLayer('Name', 'relu_3')
                      fullyConnectedLayer(64, 'Name', 'fc_4')
                      reluLayer('Name', 'relu_4')
                      fullyConnectedLayer(32, 'Name', 'fc_5')
                      reluLayer('Name', 'relu_5')
                      fullyConnectedLayer(outDim, 'Name', 'fc_6')];
            lgraph  = layerGraph(layers);
            obj.net = dlnetwork(lgraph);
            % Init params
            obj.n_states           = inDim;
            obj.n_actions          = outDim;
            obj.learn_step_counter = 0;
            obj.memory_counter     = 0;
            % Init memory
            obj.memory_scurr = zeros(inDim, obj.MEMORY_CAPACITY);
            obj.memory_a     = zeros(1, obj.MEMORY_CAPACITY);
            obj.memory_r     = zeros(1, obj.MEMORY_CAPACITY);
            obj.memory_snext = zeros(inDim, obj.MEMORY_CAPACITY);
        end
        function a = choose_action(obj, scurr)
        % Forward propagation, and choose action according to probability
            q_eval_raw = extractdata(predict(obj.net, dlarray(scurr, 'CB')));
            [~, idxs]  = sort(q_eval_raw);
            idxs       = obj.n_actions + 1 - idxs;
            % Idxs
            prob = obj.idxs2prob(idxs);
            prob = prob ./ sum(prob);
            a    = zeros(1, size(scurr, 2));
            for i = 1 : size(scurr, 2)
                a(i) = randsrc(1, 1, [1:obj.n_actions; prob(i,:)]);
            end
        end
        function prob = idxs2prob(obj, idxs)
        % Map from sorting index to probability
            prob_table = [70, 28, 10, 8];
            if idxs < 5
                prob = prob_table(idxs);
            else
                prob = 5;
            end
        end

        function obj = store_transition(obj, scurr, a, r, snext)
        % Store state transition table
            idx = rem(obj.memory_counter, obj.MEMORY_CAPACITY) + 1;
            obj.memory_scurr(:, idx) = scurr;
            obj.memory_a(:, idx)     = a;
            obj.memory_r(:, idx)     = r;
            obj.memory_snext(:, idx) = snext;
            obj.memory_counter       = obj.memory_counter + 1;
        end
        function learn(obj)
        % Train one step
            obj.learn_step_counter = obj.learn_step_counter + 1;
            idxs        = randi(min(obj.MEMORY_CAPACITY, obj.memory_counter), 1, obj.BATCH_SIZE);
            batch_scurr = obj.memory_scurr(:, idxs);
            batch_a     = obj.memory_a(:, idxs);
            batch_r     = obj.memory_r(:, idxs);
            [loss,grad] = dlfeval(@obj.modelLoss, obj.net, batch_scurr, batch_a, batch_r);
            obj.net     = dlupdate(@obj.sgdFunction, obj.net, grad);
        end
        function [loss, grad] = modelLoss(obj, net, scurr, a, r)
        % Calculate loss and gradient
            q_eval_raw = predict(net, dlarray(scurr, 'CB'));
            q_eval     = dlarray(zeros(1, obj.BATCH_SIZE), 'CB');
            for i = 1 : obj.BATCH_SIZE
                q_eval(i) = q_eval_raw(a(i), i);
            end
            loss = mse(q_eval, dlarray(r, 'CB'));
            grad = dlgradient(loss, net.Learnables);
        end
        function upd_param = sgdFunction(obj, param, grad)
        % Update parameters via SGD
            upd_param = param - obj.LR .* grad;
        end
    end
end
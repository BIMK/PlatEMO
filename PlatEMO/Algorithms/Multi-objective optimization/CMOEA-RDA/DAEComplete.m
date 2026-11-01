classdef DAEComplete < handle

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    properties
    % ----- Network structure parameters -----
    nVisible    % Visible layer dimension
    nHidden     % Hidden dimension
    
    % ----- Network weight (encoder)-----
    W1          % Encoder weight matrix [nHidden x nVisible]           
    b1          % Encoder bias [1 x nHidden]

    % ----- Network weight (decoder) -----
    W2          % Decoder weight matrix [nVisible x nHidden]
    b2          % Decoder bias [1 x nVisible]

    % ----- Data normalization parameters -----
    dataMin     % The minimum value of each dimension in the training data [1 x nVisible]
    dataMax     % The maximum value of each dimension in the training data [1 x nVisible]

    % ----- Activation function type -----
    activation  % string: 'sigmoid', 'relu', 'tanh', 'linear'
    end

    methods
        %% ====================================================================
        % Constructor: Automatically called when creating DAE objects
        % ====================================================================
        function obj = DAEComplete(nVisible, nHidden, activation)
            % input parameters:
            %   nVisible  : Input data dimension 
            %   nHidden   : Hidden dimension 
            %   activation: Activation function type 
            
            
            if nargin < 2
                nHidden = 1;  
            end

            if nargin < 3
                activation = 'sigmoid';  
            end
           
            % Save network structure parameters
            obj.nVisible   = nVisible;
            obj.nHidden    = nHidden;
            obj.activation = lower(activation);  
            
            % ----- Weight initialization (Xavier/He initialization) -----

            if strcmp(obj.activation, 'relu')
                scale = sqrt(2 / (nVisible + nHidden));
            else
                scale = sqrt(6 / (nVisible + nHidden));
            end
            
            % Encoder weight initialization
            obj.W1 = randn(nHidden, nVisible) * scale;
            obj.b1 = zeros(1, nHidden);  % Bias initialization to 0
            
            % Decoder weight initialization
            obj.W2 = randn(nVisible, nHidden) * scale;
            obj.b2 = zeros(1, nVisible);% Bias initialization to 0
        end
        
        
        %% ====================================================================
        % Training method: Train DAE using noisy data
        % ====================================================================
        function train(obj, X, epochs, lr, batchSize, noiseLevel)
            % input parameters:
            %   X          : training data [nSamples x nVisible]
            %   epochs     : training epochs 
            %   lr         : learning rate 
            %   batchSize  : batch size 
            %   noiseLevel : noise level 
            %                
            
   
            if nargin < 3, epochs = 20; end
            if nargin < 4, lr = 0.01; end
            if nargin < 5, batchSize = 20; end
            if nargin < 6, noiseLevel = 0.2; end
            
            % ----- data normalization-----
            obj.dataMin = min(X, [], 1);        
            obj.dataMax = max(X, [], 1);        
            X_norm      = (X - repmat(obj.dataMin, size(X,1), 1)) ./ ...
                          repmat(obj.dataMax - obj.dataMin, size(X,1), 1);
            
            % ----- Initialize momentum term -----
            vW1      = zeros(size(obj.W1));
            vb1      = zeros(size(obj.b1)); 
            vW2      = zeros(size(obj.W2)); 
            vb2      = zeros(size(obj.b2));  
            momentum = 0.9;  % Momentum coefficient 
           
            % ----- Main training cycle -----
            for epoch = 1 : epochs

                shuffleIdx = randperm(size(X_norm, 1));
                X_shuffled = X_norm(shuffleIdx, :);
                
                totalLoss = 0;  

                numBatches = ceil(size(X_norm, 1) / batchSize);
                
                % ----- Batch processing training -----
                for batch = 1 : numBatches
                    startIdx = (batch-1)*batchSize + 1;
                    endIdx   = min(batch*batchSize, size(X_norm, 1));
                    batchX   = X_shuffled(startIdx:endIdx, :);
                    batchSizeActual = size(batchX, 1); 
                    
                    % ----- noise addition(Masking Noise) -----
                    mask = rand(size(batchX)) > noiseLevel;
                    batchX_noisy = batchX .* mask;  
                    
                    % ----- Forward Pass -----
                    Z1 = batchX_noisy * obj.W1' + ...
                         repmat(obj.b1, batchSizeActual, 1);
                    H = obj.applyActivation(Z1); 
                    
                    % Decode
                    Z2 = H * obj.W2' + ...
                         repmat(obj.b2, batchSizeActual, 1);
                    output = obj.applyActivation(Z2);  
                    
                    % ----- calculate the loss-----
                    loss      = mean((output - batchX).^2, 'all');
                    totalLoss = totalLoss + loss * batchSizeActual;
                    
                    % ----- Backward Pass-----
                    dLoss_dOutput = 2 * (output - batchX) / batchSizeActual;
                    
                    % Decoder gradient
                    dZ2 = dLoss_dOutput .* obj.activationDerivative(output);
                    dW2 = H' * dZ2;  
                    db2 = sum(dZ2, 1); 
                    
                    % Hidden layer error 
                    dH  = dZ2 * obj.W2;  
                    dZ1 = dH .* obj.activationDerivative(H);  
                    dW1 = batchX_noisy' * dZ1; 
                    db1 = sum(dZ1, 1); 
                    
                    % ----- Gradient update with momentum-----
                    % Update decoder weights
                    vW2    = momentum * vW2 - lr * dW2';
                    vb2    = momentum * vb2 - lr * db2; 
                    vW1    = momentum * vW1 - lr * dW1';  
                    vb1    = momentum * vb1 - lr * db1;    
                    obj.W2 = obj.W2 + vW2;  
                    obj.b2 = obj.b2 + vb2;  
                    obj.W1 = obj.W1 + vW1;  
                    obj.b1 = obj.b1 + vb1;  
                end
            end
        end
       
        function H = reduce(obj, X)
            % Function: Compress input data X into low dimensional features H
            X_norm = (X - repmat(obj.dataMin, size(X,1), 1)) ./ ...
                     repmat(obj.dataMax - obj.dataMin, size(X,1), 1);
            
            % ----- Forward propagation of encoder -----
            Z = X_norm * obj.W1' + repmat(obj.b1, size(X,1), 1);
            H = obj.applyActivation(Z);  
        end

        function Y = applyActivation(obj, Z)
            % Function: Calculate the activation function value based on obj. activation
            switch obj.activation
                case 'sigmoid'
                    Y = 1 ./ (1 + exp(-Z));
                case 'relu'
                    Y = max(0, Z);
                case 'tanh'
                    Y = tanh(Z);
                case 'linear'
                    Y = Z;
                otherwise
                    error('Unknown activation function: %s', obj.activation);
            end
        end
        
        function dY = activationDerivative(obj, Y)
            % Function: Calculate the derivative of the activation function (for backpropagation)
            switch obj.activation
                case 'sigmoid'
                    dY = Y .* (1 - Y);
                case 'relu'
                    dY = double(Y > 0);
                case 'tanh'
                    dY = 1 - Y.^2;
                case 'linear'
                    dY = ones(size(Y));
                otherwise
                    error('Unknown activation function: %s', obj.activation);
            end
        end
        
    end 
end 
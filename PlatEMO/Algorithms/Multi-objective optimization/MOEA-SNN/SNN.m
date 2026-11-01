classdef SNN < handle

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: mingfei@westlake.edu.cn)

    properties
        nInput;
        nHidden;
        Epoch;
        LearnRate;
        Weights;
        Bias;
        Decay;
        Threshold;
    end
    methods
        %% Constructor
        function obj = SNN(nInput, nHidden, Epoch, LearnRate)
            obj.nInput    = nInput;
            obj.nHidden   = nHidden;
            obj.Epoch     = Epoch;
            obj.LearnRate = LearnRate;
            obj.Weights   = 0.1*randn(nInput,nHidden,'single');
            obj.Bias      = zeros(1,nHidden,'single');
            obj.Decay     = 0.95;
            obj.Threshold = 0.5;
        end
        %% Train
        function train(obj,X)
            N          = size(X,1);
            T          = 10;                           
            targetRate = obj.nHidden/obj.nInput;                
            lambda_bp  = 0.1;                  
            lambda_s   = 0.2;                    
            rho        = 0.05;                       
            epsVal     = 1e-8;
            if isscalar(obj.Threshold)
                obj.Threshold = obj.Threshold * ones(1,obj.nHidden,'single');
            end
            for epoch = 1 : obj.Epoch
                idx = randperm(N);
                X   = X(idx,:);
                V   = zeros(N, obj.nHidden, 'single');

                for t = 1 : T
                    I = X * obj.Weights + obj.Bias;

                    V = obj.Decay * V + I;

                    spike    = V > obj.Threshold;
                    V(spike) = 0;

                    obj.Weights = obj.Weights + obj.LearnRate * (X' * single(spike)) / N;
                    obj.Bias    = obj.Bias + obj.LearnRate * mean(single(spike),1);
                end

                H    = obj.sigmoid(X * obj.Weights + obj.Bias);  
                Xhat = obj.sigmoid(H * obj.Weights');         

                kl     = X .* log((X + epsVal) ./ (Xhat + epsVal)) + ...
                         (1 - X) .* log((1 - X + epsVal) ./ (1 - Xhat + epsVal));
                klLoss = mean(sum(kl,2));

                rho_hat    = mean(H,1);   
                sparseLoss = rho .* log((rho + epsVal) ./ (rho_hat + epsVal)) + ...
                             (1 - rho) .* log((1 - rho + epsVal) ./ (1 - rho_hat + epsVal));
                sparseLoss = sum(sparseLoss);

                dXhat = (Xhat - X);            
                dH    = dXhat * obj.Weights .* H .* (1 - H);  

                dSparse = lambda_s * ((-rho ./ (rho_hat + epsVal)) + ...
                          (1 - rho) ./ (1 - rho_hat + epsVal)) / N;
                dH      = dH + ones(N,1) * dSparse;                     

                gradW = (X' * dH + (H' * dXhat)') / N;
                gradB = mean(dH,1);

                obj.Weights = obj.Weights - lambda_bp * obj.LearnRate * gradW;
                obj.Bias    = obj.Bias    - lambda_bp * obj.LearnRate * gradB;

                obj.Weights = max(min(obj.Weights,1.0),-1.0);

                neuronRate    = mean(spike(:),1);
                obj.Threshold = obj.Threshold + 0.1 * (neuronRate - targetRate);
            end
        end
        %% Reduce
        function H = reduce(obj,X)
            V = X * obj.Weights + obj.Bias;
            H = V > obj.Threshold; 
        end
        %% Recover
        function Xhat = recover(obj,H)
            V    = H * obj.Weights' + zeros(1,obj.nInput,'like',H);
            thr  = mean(obj.Threshold(:));
            Xhat = V > thr; 
        end
        %% Sigmoid
        function y = sigmoid(~, x)
            y = 1 ./ (1 + exp(-x));
        end
    end
end
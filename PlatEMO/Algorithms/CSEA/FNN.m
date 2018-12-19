classdef FNN < handle
% Feedforward neural network

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        nHidden = 0;    % Size of the hidden layer
        nIter   = 0;    % Number of iterations for each learning
        WA      = [];	% Weights between the input layper and the hidden layer
        WB      = [];	% Weights between the hidden layer and the output layer
    end
    methods
        %% Constructor
        function obj = FNN(nHidden,nIter)
            obj.nHidden = nHidden;
            obj.nIter   = nIter;
        end
        
        %% Train
        function train(obj,X,T)
            if isempty(obj.WA) || isempty(obj.WB)
                obj.WA = randn(size(X,2)+1,obj.nHidden)./100;
                obj.WB = randn(obj.nHidden+1,size(T,2))./100;
            end
            miu = 1;
            k   = 2;
            for iter = 1 : obj.nIter
                % Calculate the predictive output
                [Z,Y] = obj.predict(X);
                MSE   = mse(Z-T);
                % Calculate the Jacobian matrix
                J = zeros(size(X,1),numel(obj.WA)+numel(obj.WB));
                for i = 1 : size(X,1)
                    P      = Z(i,:).*(1-Z(i,:));
                    Q      = P*obj.WB(2:end,:)'.*Y(i,:).*(1-Y(i,:));
                    DB     = [1,Y(i,:)]'*P;
                    DA     = [1,X(i,:)]'*Q;
                    J(i,:) = [DA(:);DB(:)];
                end
                % Update the value of each weight
                J = roundn(J,-4);
                while true
                    Delta = -(J'*J+miu*eye(size(J,2)))^-1*J'*sum(Z-T,2);
                    newWA = obj.WA + reshape(Delta(1:numel(obj.WA)),size(obj.WA));
                    newWB = obj.WB + reshape(Delta(numel(obj.WA)+1:end),size(obj.WB));
                    newY  = 1./(1+exp(-[ones(size(X,1),1),X]*newWA));
                    newZ  = 1./(1+exp(-[ones(size(Y,1),1),newY]*newWB));
                    if MSE < 1e-4
                        return;
                    end
                    if mse(newZ-T) < MSE
                        obj.WA = newWA;
                        obj.WB = newWB;
                        miu    = miu/k;
                        break;
                    elseif miu > 1e4
                        return;
                    else
                        miu = miu*k;
                    end
                end
            end
        end
        
        %% Predict
        function [Z,Y] = predict(obj,X)
            Y = 1./(1+exp(-[ones(size(X,1),1),X]*obj.WA));
            Z = 1./(1+exp(-[ones(size(Y,1),1),Y]*obj.WB));
        end
    end
end
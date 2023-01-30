classdef DAE < handle
% Feedforward neural network

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        nVisible  = 0;
        nHidden   = 0;
        Epoch     = 10;
        BatchSize = 1;
        InputZeroMaskedFraction = 0.5;
        Momentum  = 0.5;
        LearnRate = 0.1;
        WA        = [];
        WB        = [];
        lower     = [];
        upper     = [];
    end
    methods
        %% Constructor
        function obj = DAE(nVisible,nHidden,Epoch,BatchSize,InputZeroMaskedFraction,Momentum,LearnRate)
            obj.nVisible  = nVisible;
            obj.nHidden   = nHidden;
            obj.Epoch     = Epoch;
            obj.BatchSize = BatchSize;
            obj.InputZeroMaskedFraction = InputZeroMaskedFraction;
            obj.Momentum  = Momentum;
            obj.LearnRate = LearnRate; 
            obj.WA = (rand(nHidden,nVisible+1)-0.5)*8*sqrt(6/(nHidden+nVisible));   
            obj.WB = (rand(nVisible,nHidden+1)-0.5)*8*sqrt(6/(nVisible+nHidden));  
        end
        %% Train
        function train(obj,X)
            obj.lower = min(X,[],1);
            obj.upper = max(X,[],1);
            X = (X-repmat(obj.lower,size(X,1),1))./repmat(obj.upper-obj.lower,size(X,1),1);
            vW{1} = zeros(size(obj.WA));
            vW{2} = zeros(size(obj.WB));
            if(obj.InputZeroMaskedFraction ~= 0)
                theta = rand(size(X)) > obj.InputZeroMaskedFraction;
            else
                theta = true(size(X));
            end
            X_temp = X.*theta;
            X_temp = [ones(size(X,1),1),X_temp];
            for i = 1 : obj.Epoch
                kk = randperm(size(X,1));
                for batch = 1 : size(X,1)/obj.BatchSize
                    batch_x = X_temp(kk((batch-1)*obj.BatchSize+1:batch*obj.BatchSize),:);
                    batch_y = X(kk((batch-1)*obj.BatchSize+1:batch*obj.BatchSize),:);

                    % Feedforward pass
                    poshid1 = 1./(1+exp(-batch_x*obj.WA'));
                    poshid1 = [ones(obj.BatchSize,1),poshid1];
                    poshid2 = 1./(1+exp(-poshid1*obj.WB'));

                    % BP
                    e     = batch_y - poshid2;
                    d{3}  = -e.*(poshid2.*(1-poshid2));
                    d_act = poshid1.*(1-poshid1);
                    d{2}  = d{3}*obj.WB.*d_act;
                    for i = 1 : 2
                        if i+1 == 3
                            dW{i} = (d{i+1}'*poshid1/size(d{3},1));
                        else
                            dW{i} = (d{i+1}(:,2:end)'*batch_x)/size(d{i+1},1);
                        end
                    end
                    for i = 1 : 2
                        dW{i} = obj.LearnRate*dW{i};
                        if obj.Momentum > 0
                            vW{i} = obj.Momentum*vW{i} + dW{i};
                            dW{i} = vW{i};
                        end
                        if i == 1
                            obj.WA = obj.WA - dW{i};
                        else
                            obj.WB = obj.WB - dW{i};
                        end
                    end
                end
            end
        end
        %% Reduce
        function H = reduce(obj,X)
            X = (X-repmat(obj.lower,size(X,1),1))./repmat(obj.upper-obj.lower,size(X,1),1);
            H = 1./(1+exp(-X*obj.WA(:,2:end)'-repmat(obj.WA(:,1)',size(X,1),1)));
        end
        %% Recover
        function X = recover(obj,H)
            X = 1./(1+exp(-H*obj.WB(:,2:end)'-repmat(obj.WB(:,1)',size(H,1),1)));
            X = X.*repmat(obj.upper-obj.lower,size(X,1),1) + repmat(obj.lower,size(X,1),1);
        end
    end
end
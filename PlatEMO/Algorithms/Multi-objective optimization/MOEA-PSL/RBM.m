classdef RBM < handle
% Restricted Boltzmann Machine

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
        Penalty   = 0.01;
        Momentum  = 0.5;
        LearnRate = 0.1;
        Weight    = [];
        vBias     = [];
        hBias     = [];
    end
    methods
        %% Constructor
        function obj = RBM(nVisible,nHidden,Epoch,BatchSize,Penalty,Momentum,LearnRate)
            obj.nVisible  = nVisible;
            obj.nHidden   = nHidden;
            obj.Epoch     = Epoch;
            obj.BatchSize = BatchSize;
            obj.Penalty   = Penalty;
            obj.Momentum  = Momentum;
            obj.LearnRate = LearnRate;
            obj.Weight = 0.1 * randn(obj.nVisible,obj.nHidden);
            obj.vBias  = zeros(1,obj.nVisible);
            obj.hBias  = zeros(1,obj.nHidden);
        end
        %% Train
        function train(obj,X)
            vishidinc  = zeros(size(obj.Weight));
	        hidbiasinc = zeros(size(obj.hBias));
	        visbiasinc = zeros(size(obj.vBias));
            for epoch = 1 : obj.Epoch
                if obj.Epoch > 5
                    obj.Momentum = 0.9;
                end
                kk = randperm(size(X,1));
                for batch = 1 : size(X,1)/obj.BatchSize
                    batchdata = X(kk((batch-1)*obj.BatchSize+1:batch*obj.BatchSize),:);

                    % Positive phase
                    poshidprobs  = 1./(1+exp(-batchdata*obj.Weight-repmat(obj.hBias,obj.BatchSize,1))); 
                    poshidstates = poshidprobs > rand(obj.BatchSize,obj.nHidden);

                    % Negative phase
                    negdataprobs = 1./(1+exp(-poshidstates*obj.Weight'-repmat(obj.vBias,obj.BatchSize,1)));
                    negdata      = negdataprobs > rand(obj.BatchSize,obj.nVisible);
                    neghidprobs  = 1./(1+exp(-negdata*obj.Weight-repmat(obj.hBias,obj.BatchSize,1))); 

                    % Update weight
                    posprods   = batchdata' * poshidprobs;
                    negprods   = negdataprobs' * neghidprobs;
                    poshidact  = sum(poshidprobs);
		            posvisact  = sum(batchdata);
                    neghidact  = sum(neghidprobs);
		            negvisact  = sum(negdata); 
                    vishidinc  = obj.Momentum*vishidinc + obj.LearnRate*(((posprods-negprods)/obj.BatchSize)-obj.Penalty*obj.Weight);
                    visbiasinc = obj.Momentum*visbiasinc + (obj.LearnRate/obj.BatchSize)*(posvisact-negvisact);
                    hidbiasinc = obj.Momentum*hidbiasinc + (obj.LearnRate/obj.BatchSize)*(poshidact-neghidact);
                    obj.Weight = obj.Weight + vishidinc;
                    obj.vBias  = obj.vBias + visbiasinc;
                    obj.hBias  = obj.hBias + hidbiasinc;
                end                
            end
        end
        %% Reduce
        function H = reduce(obj,X)
            H = 1./(1+exp(-X*obj.Weight-repmat(obj.hBias,size(X,1),1))) > rand(size(X,1),size(obj.Weight,2));
        end
        %% Recover
        function X = recover(obj,H)
            X = 1./(1+exp(-H*obj.Weight'-repmat(obj.vBias,size(H,1),1))) > rand(size(H,1),size(obj.Weight,1));
        end
    end
end
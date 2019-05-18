classdef Sparse_NN < PROBLEM
% <problem> <Sparse MOP>
% The neural network training problem (only for bi-category classification)
% dataNo  ---  1 --- Number of dataset
% nHidden --- 20 --- Size of hidden layer

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zhang, C. Wang, and Y. Jin, An evolutionary algorithm for
% large-scale sparse multi-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% The datasets are taken from the UCI machine learning repository in
% http://archive.ics.uci.edu/ml/index.php
% No.   Name                              Samples Features Classes
% 1     Wine                                178      13       3
% 2     Statlog_Australian                  690      14       2
% 3     Climate                             540      18       2
% 4     Parkinsons                          195      22       2
% 5     Statlog_German                     1000      24       2
% 6     Breast_cancer_Wisconsin_Diagnostic  569      30       2
% 7     Ionosphere                          351      34       2
% 8     SPECTF_Heart                        267      44       2
% 9     Lung_cancer                          32      56       3
% 10    Connectionist_bench_Sonar           208      60       2
% 11    Libras_movement                     360      90      15
% 12    Hill_Valley                         606     100       2
% 13    MUSK1                               476     166       2
% 14    Semeion_handwritten_digit          1593     256      10
% 15    LSVT_voice_rehabilitation           126     310       2
% 16    Madelon                            2600     500       2
% 17    ISOLET                             1557     617      26
% 18    Multiple_features                  2000     649      10
% 19    CNAE9                              1080     857       9

    properties(Access = private)
        TrainIn;    % Input of training set
        TrainOut;   % Output of training set
        TestIn;   	% Input of validation set
        TestOut;  	% Output of validation set
        nHidden;    % Size of hidden layer
    end
    methods
        function obj = Sparse_NN()
            % Load data
            [dataNo,obj.nHidden] = obj.Global.ParameterSet(1,20);
            str    = {'Wine','Statlog_Australian','Climate','Parkinsons','Statlog_German','Breast_cancer_Wisconsin_Diagnostic',...
                      'Ionosphere','SPECTF_Heart','Lung_cancer','Connectionist_bench_Sonar','Libras_movement','Hill_Valley',...
                      'MUSK1','Semeion_handwritten_digit','LSVT_voice_rehabilitation','Madelon','ISOLET','Multiple_features','CNAE9'};
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_FS_NN.mat'),'Dataset');
            Data = Dataset.(str{dataNo});
            Mean = mean(Data(:,1:end-1),1);
            Std  = std(Data(:,1:end-1),[],1);
            Data(:,1:end-1) = (Data(:,1:end-1)-repmat(Mean,size(Data,1),1))./repmat(Std,size(Data,1),1);
            Data(:,end)     = Data(:,end) == Data(1,end);	% Only for bi-category classification
            obj.TrainIn     = Data(1:ceil(end*0.8),1:end-1);
            obj.TrainOut    = Data(1:ceil(end*0.8),end);
            obj.TestIn      = Data(ceil(end*0.8)+1:end,1:end-1);
            obj.TestOut     = Data(ceil(end*0.8)+1:end,end);
            % Parameter setting
            obj.Global.M        = 2;
            obj.Global.D        = (size(obj.TrainIn,2)+1)*obj.nHidden + (obj.nHidden+1)*size(obj.TrainOut,2);
            obj.Global.lower    = zeros(1,obj.Global.D) - 1;
            obj.Global.upper    = zeros(1,obj.Global.D) + 1;
            obj.Global.encoding = 'real';
        end
        %% Generate initial population
        function PopDec = Init(obj,N)
            PopDec = (rand(N,obj.Global.D)-0.5)*2.*randi([0 1],N,obj.Global.D);
        end
        %% Fine-tune solutions
        function PopDec = CalDec(obj,PopDec)
            for i = 1 : size(PopDec,1)
                W1          = reshape(PopDec(i,1:(size(obj.TrainIn,2)+1)*obj.nHidden),size(obj.TrainIn,2)+1,obj.nHidden);
                W2          = reshape(PopDec(i,(size(obj.TrainIn,2)+1)*obj.nHidden+1:end),obj.nHidden+1,size(obj.TrainOut,2));
                [W1,W2]     = Train(obj.TrainIn,obj.TrainOut,W1,W2,1);
                PopDec(i,:) = [W1(:)',W2(:)'];
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),2);
            for i = 1 : size(PopDec,1)
                W1 = reshape(PopDec(i,1:(size(obj.TrainIn,2)+1)*obj.nHidden),size(obj.TrainIn,2)+1,obj.nHidden);
                W2 = reshape(PopDec(i,(size(obj.TrainIn,2)+1)*obj.nHidden+1:end),obj.nHidden+1,size(obj.TrainOut,2));
                Z  = Predict(obj.TrainIn,W1,W2);
                PopObj(i,1) = mean(PopDec(i,:)~=0);
                PopObj(i,2) = mean(round(Z)~=obj.TrainOut);
            end
        end
        %% Draw special figure
        function Draw(obj,PopDec)
            PopObj = zeros(size(PopDec,1),2);
            for i = 1 : size(PopDec,1)
                W1 = reshape(PopDec(i,1:(size(obj.TestIn,2)+1)*obj.nHidden),size(obj.TestIn,2)+1,obj.nHidden);
                W2 = reshape(PopDec(i,(size(obj.TestIn,2)+1)*obj.nHidden+1:end),obj.nHidden+1,size(obj.TestOut,2));
                Z  = Predict(obj.TestIn,W1,W2);
                PopObj(i,1) = mean(PopDec(i,:)~=0);
                PopObj(i,2) = mean(round(Z)~=obj.TestOut);
            end
            cla; Draw(PopObj);
            xlabel('Complexity'); ylabel('Test error');
        end
    end
end

function [W1,W2] = Train(X,T,W1,W2,nEpoch)
    for epoch = 1 : nEpoch
        [Z,Y] = Predict(X,W1,W2);
        P     = (Z-T).*Z.*(1-Z);
        Q     = P*W2(2:end,:)'.*(1-Y.^2);
        D1    = 0;
        D2    = 0;
        for i = 1 : size(X,1)
            D2 = D2 + [0,Y(i,:)]'*P(i,:);
            D1 = D1 + [0,X(i,:)]'*Q(i,:);
        end
        W1 = W1 - D1/size(X,1);
        W2 = W2 - D2/size(X,1);
    end
end

function [Z,Y] = Predict(X,W1,W2)
    Y = 1 - 2./(1+exp(2*[zeros(size(X,1),1),X]*W1));
    Z = 1./(1+exp(-[zeros(size(Y,1),1),Y]*W2));
end
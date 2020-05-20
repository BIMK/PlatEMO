classdef Sparse_IS < PROBLEM
% <problem> <Sparse MOP>
% The instance selection problem
% dataNo --- 1 --- Number of dataset

%------------------------------- Reference --------------------------------
% Y. Tian, C. Lu, X. Zhang, K. C. Tan, and Y. Jin, Solving large-scale
% multi-objective optimization problems with sparse optimal solutions via
% unsupervised neural networks, IEEE Transactions on Cybernetics, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% The datasets are taken from the LIBSVM Data in
% https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html
% and the UCI machine learning repository in
% http://archive.ics.uci.edu/ml/index.php
% No.   Name      	Samples Features Classes
% 1     Fourclass	  862       3       2
% 2     Abalone  	 4177       9       2
% 3     Phishing	11055      68       2

    properties(Access = private)
        Data;   % Dataset
    end
    methods
        %% Initialization
        function obj = Sparse_IS()
            % Load data
            dataNo    = obj.Global.ParameterSet(1);
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_IS.mat'),'Dataset');
            str       = {'Fourclass','Abalone','Phishing'};
            obj.Data  = Dataset.(str{dataNo});
            if issparse(obj.Data)
                obj.Data = full(obj.Data);
            end
            [Temp,~] = mapminmax(obj.Data(:,1:end-1)', 0, 1);
            obj.Data(:,1:end-1) = Temp';
            % Parameter setting
            obj.Global.M = 2;
            obj.Global.D = size(obj.Data,1);
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'binary';
        end
		
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec = logical(round(PopDec));          
            PopObj = zeros(size(PopDec,1),obj.Global.M);
            for i = 1 : size(PopObj,1)
                TrainIn  = obj.Data(PopDec(i,:),1:end-1);
                TrainOut = obj.Data(PopDec(i,:),end);
                ValidIn  = obj.Data(~PopDec(i,:),1:end-1);
                ValidOut = obj.Data(~PopDec(i,:),end);
                if mean(PopDec(i,:)) == 1
                    PopObj(i,1) = 1;
                    PopObj(i,2) = 0;
                elseif length(unique(TrainOut)) == 1
                    PopObj(i,1) = mean(PopDec(i,:));
                    PopObj(i,2) = 1;
                elseif mean(PopDec(i,:)) > 0
                    SVMModel = fitcsvm(TrainIn,TrainOut);
                    label = predict(SVMModel,ValidIn);
                    PopObj(i,1) = mean(PopDec(i,:));
                    PopObj(i,2) = mean(label~=ValidOut);
                else
                    PopObj(i,1) = 0;
                    PopObj(i,2) = 1;
                end
            end
        end
    end
end
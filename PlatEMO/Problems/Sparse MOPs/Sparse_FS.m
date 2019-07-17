classdef Sparse_FS < PROBLEM
% <problem> <Sparse MOP>
% The feature selection problem
% dataNo --- 1 --- Number of dataset

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
% 1     MUSK1                               476     166       2
% 2     Semeion_handwritten_digit          1593     256      10
% 3     LSVT_voice_rehabilitation           126     310       2
% 4     ISOLET                             1557     617      26

    properties(Access = private)
        TrainIn;    % Input of training set
        TrainOut;   % Output of training set
        ValidIn;    % Input of validation set
        ValidOut;   % Output of validation set
        Category;   % Output label set
    end
    methods
        %% Initialization
        function obj = Sparse_FS()
            % Load data
            dataNo = obj.Global.ParameterSet(1);
            str    = {'MUSK1','Semeion_handwritten_digit','LSVT_voice_rehabilitation','ISOLET'};
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_FS.mat'),'Dataset');
            Data = Dataset.(str{dataNo});
            Fmin = min(Data(:,1:end-1),[],1);
            Fmax = max(Data(:,1:end-1),[],1);
            Data(:,1:end-1) = (Data(:,1:end-1)-repmat(Fmin,size(Data,1),1))./repmat(Fmax-Fmin,size(Data,1),1);
            obj.Category    = unique(Data(:,end));
            obj.TrainIn     = Data(1:ceil(end*0.8),1:end-1);
            obj.TrainOut    = Data(1:ceil(end*0.8),end);
            obj.ValidIn     = Data(ceil(end*0.8)+1:end,1:end-1);
            obj.ValidOut    = Data(ceil(end*0.8)+1:end,end);
            % Parameter setting
            obj.Global.M        = 2;
            obj.Global.D        = size(obj.TrainIn,2);
            obj.Global.encoding = 'binary';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec   = logical(PopDec);
            PopObj   = zeros(size(PopDec,1),2);
            for i = 1 : size(PopObj,1)
                [~,Rank] = sort(pdist2(obj.ValidIn(:,PopDec(i,:)),obj.TrainIn(:,PopDec(i,:))),2);
                [~,Out]  = max(hist(obj.TrainOut(Rank(:,1:3))',obj.Category),[],1);
                Out      = obj.Category(Out);
            	PopObj(i,1) = mean(PopDec(i,:));
                PopObj(i,2) = mean(Out~=obj.ValidOut);
            end
        end
    end
end
classdef TREE4 < PROBLEM
% <multi> <real> <large/none> <constrained> <expensive/none>
% The time-varying ratio error estimation problem
% T --- 1000 --- Length of data (related to the number of variables)

%------------------------------- Reference --------------------------------
% C. He, R. Cheng, C. Zhang, Y. Tian, Q. Chen, and X. Yao, Evolutionary
% large-scale multiobjective optimization for ratio error estimation of
% voltage transformers, IEEE Transactions on Evolutionary Computation,
% 2020, 24(5): 868-881.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Data;   % Dataset
        Mean;   % Mean values of the dataset
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Load data
            T = obj.ParameterSet(1000);
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_TREE.mat'),'Dataset');
            obj.Data = Dataset.TREE4(1:min(max(T,10),end),:);
            % Set the numbers of objectives and decision variables
            K = 6;
            obj.M = 2;
            obj.D = T*K;
            % Set the upper and lower boundaries
            Lower    = zeros(T,K);
            Upper    = zeros(T,K);
            obj.Mean = zeros(T,K);
            for i = 1 : 3
                Lower(:,i)      = min(obj.Data(:,i:3:end*0.5),[],2)/1.01;
                Upper(:,i)      = max(obj.Data(:,i:3:end*0.5),[],2)/0.99;
                Lower(:,i+3)    = min(obj.Data(:,end*0.5+i:3:end),[],2)/1.01;
                Upper(:,i+3)    = max(obj.Data(:,end*0.5+i:3:end),[],2)/0.99;
                obj.Mean(:,i)   = mean(obj.Data(:,i:3:end*0.5),2);
                obj.Mean(:,i+3) = mean(obj.Data(:,end*0.5+i:3:end),2);
            end
            % The decision variables are the offset of the mean values of
            % the dataset
            obj.lower    = reshape(Lower-obj.Mean,1,[]);
            obj.upper    = reshape(Upper-obj.Mean,1,[]);
            obj.encoding = ones(1,obj.D);
        end
        %% Generate initial solutions
        function Population = Initialization(obj,N)
            if nargin < 2; N = obj.N; end
            PopDec = obj.Mean(:)'.*(rand(N,obj.D)*0.008-0.004);
            PopDec = min(max(PopDec,repmat(obj.lower,N,1)),repmat(obj.upper,N,1));
            Population = obj.Evaluation(PopDec);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            N      = size(PopDec,1);
            KP     = size(obj.Data,2);
            [T,K]  = size(obj.Mean);
            PopDec = reshape(PopDec,N,T,K) + repmat(reshape(obj.Mean,1,T,K),N,1,1);
            PopObj = zeros(N,obj.M);
            % First objective
            eA1 = abs(repmat(reshape(obj.Data(:,1:end*0.5),1,T,KP/2),N,1,1)./repmat(PopDec(:,:,1:3),1,1,KP/6)-1);
            eA2 = abs(repmat(reshape(obj.Data(:,end*0.5+1:end),1,T,KP/2),N,1,1)./repmat(PopDec(:,:,4:6),1,1,KP/6)-1);
            eA  = cat(3,eA1,eA2);
            PopObj(:,1) = sum(sum(eA,3),2);
            % Second objective
            Delta = std(eA(:,2:end,:)-eA(:,1:end-1,:),0,2);
            PopObj(:,2) = sum(reshape(Delta,N,[]),2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            N      = size(PopDec,1);
            [T,K]  = size(obj.Mean);
            PopDec = reshape(PopDec,N,T,K) + repmat(reshape(obj.Mean,1,T,K),N,1,1);
            PopCon = TREE_CalCon(PopDec,2);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            X = zeros(1,obj.D);
            X(1:2:end) = obj.lower(1:2:end);
            X(2:2:end) = obj.upper(2:2:end);
            R = obj.CalObj(X);
        end
    end
end
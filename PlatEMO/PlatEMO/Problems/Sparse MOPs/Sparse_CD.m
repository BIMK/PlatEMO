classdef Sparse_CD < PROBLEM
% <problem> <Sparse MOP>
% The community detection problem
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

% The datasets are taken from the Alex Arenas Website in
% http://deim.urv.cat/~alexandre.arenas/data/welcome.htm
% No.   Name	Nodes   Edges
% 1     Email	1133    5451

    properties(Access = private)
        Data;	% Dataset
        MaxKKM;	% Maximum and minimum objective values for normalization
        MinKKM;
        MaxRC;
        MinRC;
    end
    methods
    	%% Initialization
        function obj = Sparse_CD()
            % Load data
            dataNo    = obj.Global.ParameterSet(1);
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_CD.mat'),'Dataset');
            str       = {'Email'};
            obj.Data  = Dataset.(str{dataNo});
            % Parameter setting
            obj.Global.M = 2;
            obj.Global.D = size(obj.Data.Adj,2);
            obj.Global.lower    = zeros(1, obj.Global.D);
            obj.Global.upper    = ones(1, obj.Global.D);
            obj.Global.encoding = 'binary';
            
            temp1 = zeros(1,obj.Global.D);
            temp1(1) = 1;
            temp2 = ones(1,obj.Global.D);           
            [obj.MaxKKM,obj.MinRC] = FunValue(temp1,obj.Data.Adj,obj.Data.ACT,obj.Global.D);
            [obj.MinKKM,obj.MaxRC] = FunValue(temp2,obj.Data.Adj,obj.Data.ACT,obj.Global.D);  
        end
		
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec = logical(round(PopDec));        
            PopObj = zeros(size(PopDec,1),obj.Global.M);
            for i = 1 : size(PopObj,1)             
                [PopObj(i,1),PopObj(i,2)]= FunValue(PopDec(i,:),obj.Data.Adj,obj.Data.ACT,obj.Global.D);                
            end
            % Normalization
            PopObj(:,1) = (PopObj(:,1)-repmat(obj.MinKKM,size(PopObj,1),1))./(repmat(obj.MaxKKM,size(PopObj,1),1)-repmat(obj.MinKKM,size(PopObj,1),1));
            PopObj(:,2) = (PopObj(:,2)-repmat(obj.MinRC,size(PopObj,1),1))./(repmat(obj.MaxRC,size(PopObj,1),1)-repmat(obj.MinRC,size(PopObj,1),1));
        end
    end
end

function  [KKM,RC] = FunValue(PopDec,Adj,ACT,M)
    %% Decoding
    PopDec      = logical(PopDec);
    Community   = 1 : M;
    IndexMatrix = 1 : M;
    if sum(PopDec) >= 2        
        [U,IndexMatrix1] = Myfunction_random(Community(1,PopDec),PopDec,IndexMatrix,ACT);    
        U = U';
        IndexMatrix1 = IndexMatrix1';
        [U_max,U_index] = max(U,[],2);      
        Utemp = U_max;     
        U_partition = U >= Utemp;  
        error_index = find(sum(U_partition,2)>1);
        U_partition(error_index,:) = 0;    
        for lp = 1 : length(error_index)
            U_partition(error_index(lp),U_index(error_index(lp))) = 1;        
        end    
        C_Community = num2cell(Community(1,PopDec));    
        for jj = 1 : size(U,2)
             C_Community{jj} = [C_Community{jj} IndexMatrix1(logical(U_partition(:,jj)),:)'];
        end 
    else
        C_Community = {IndexMatrix};
    end    
    
    %% Objectives: KKM and RC
    KKM = compute_objective(Adj,C_Community);   
    RC  = community_score(Adj,C_Community);
end

function [U_new,IndexMatrix] = Myfunction_random(Community,PopDec,IndexMatrix,ACT)
    [~,n] = size(Community);
    tmp   = ACT(PopDec,~PopDec);
    U_new = tmp./(ones(n,1)*sum(tmp));
    IndexMatrix = IndexMatrix(~PopDec);
end

function cf = compute_objective(Adj, clu_assignment)
% Compute the community fittness for each partition.
% Adj: the adjacency matrix of the network.
% clu_assignment: the cluster label vector.

    cf = 0;
    ec = 0;
    numVar  = size(Adj,1);
    clu_num = length(clu_assignment);
    for i = 1 : clu_num
        s_index = clu_assignment{i};
        s = Adj(s_index,s_index);
        s_cardinality = length(s_index);
        if s_cardinality > 0
            kins_sum = 0;
            for j = 1:s_cardinality
                kins = sum(s(j,:));
                kins_sum = kins_sum + kins;
                ec = kins_sum;
            end
            cf_s = ec*1.0/(s_cardinality);
            cf   = cf + cf_s;
        end
    end
    cf = 2*(numVar-length(clu_assignment)) - cf;
end

function cs = community_score(Adj,clu_assignment)
% Compute the community score for each partition.
% Adj: the adjacency matrix of the network.
% clu_assignment: the cluster label vector.

    cs = 0;
    de = 0;
    clu_num = length(clu_assignment);
    for i = 1 : clu_num
        s_index =clu_assignment{i};
        s = Adj(s_index,s_index);
        s_cardinality = length(s_index);
        if s_cardinality > 0
            kins_sum  = 0;
            kouts_sum = 0;
            for j = 1:s_cardinality
                kins      = sum(s(j,:));
                ksum      = sum(Adj(s_index(j),:));
                kouts     = ksum - kins;
                kins_sum  = kins_sum + kins;
                kouts_sum = kouts_sum + kouts;
                de = kouts_sum;
            end
            cf_s = de*1.0/(s_cardinality);
            cs   = cs + cf_s;
        end
    end
end
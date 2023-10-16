classdef Sparse_CD < PROBLEM
% <multi> <binary> <large/none> <expensive/none> <sparse/none>
% The community detection problem
% dataNo --- 1 --- Number of dataset

%------------------------------- Reference --------------------------------
% Y. Tian, C. Lu, X. Zhang, K. C. Tan, and Y. Jin, Solving large-scale
% multi-objective optimization problems with sparse optimal solutions via
% unsupervised neural networks, IEEE Transactions on Cybernetics, 2021,
% 51(6): 3115-3128.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% The datasets are taken from
% http://www-personal.umich.edu/~mejn/netdata/
% No.   Name        Nodes   Edges
% 1     Karate      34      78
% 2     Dolphin     62      159
% 3     Polbook     105     441
% 4     Football    115     613

    properties(Access = private)
        Adj;    % Adjacency matrix of the network
        ACT;    % Similarity between nodes
        G;      % The graph object
        MaxKKM;	% Maximum and minimum objective values for normalization
        MinKKM;
        MaxRC;
        MinRC;
    end
    methods
    	%% Default settings of the problem
        function Setting(obj)
            % Load data
            dataNo    = obj.ParameterSet(1);
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_CD.mat'),'Dataset');
            str     = {'Karate','Dolphin','Polbook','Football'};
            obj.Adj = Dataset.(str{dataNo});
            obj.ACT = random_walk_distance(obj.Adj).^(-2);
            obj.G   = graph(obj.Adj);
            % Parameter setting
            obj.M = 2;
            obj.D = size(obj.Adj,2);
            obj.lower    = zeros(1, obj.D);
            obj.upper    = ones(1, obj.D);
            obj.encoding = 4 + zeros(1,obj.D);
            % Maximum and minimum objective values for normalization
            C = Decoding([1,zeros(1,obj.D-1)],obj.ACT,obj.D);
            obj.MaxKKM = KKM(obj.Adj,C);
            obj.MinRC  = RC(obj.Adj,C);
            C = Decoding(ones(1,obj.D),obj.ACT,obj.D);
            obj.MinKKM = KKM(obj.Adj,C);
            obj.MaxRC  = RC(obj.Adj,C);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec = logical(round(PopDec));        
            PopObj = zeros(size(PopDec,1),obj.M);
            for i = 1 : size(PopObj,1)
                C = Decoding(PopDec(i,:),obj.ACT,obj.D);
                PopObj(i,1) = KKM(obj.Adj,C);
                PopObj(i,2) = RC(obj.Adj,C);
            end
            PopObj(:,1) = (PopObj(:,1)-obj.MinKKM)./(obj.MaxKKM-obj.MinKKM);
            PopObj(:,2) = (PopObj(:,2)-obj.MinRC)./(obj.MaxRC-obj.MinRC);
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            C = cell(1,length(Population));
            Q = zeros(1,length(Population));
            for i = 1 : length(Population)
                C{i} = Decoding(Population(i).dec,obj.ACT,obj.D);
                Q(i) = Modularity(obj.Adj,C{i});
            end
            [~,best] = max(Q);
            h = plot(Draw([]),obj.G,'MarkerSize',6,'EdgeColor',[.3 .3 .3]);
            tempStream = RandStream('mlfg6331_64','Seed',2);
            for i = 1 : length(C{best})
                highlight(h,C{best}{i},'NodeColor',rand(tempStream,1,3));
            end
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            Draw(Population.objs,{'Kernel k-means','Ratio cut',[]});
        end
    end
end

function ACT = random_walk_distance(Adj)
% Random walk based distance

	n       = length(Adj);
    D       = diag(sum(Adj,2));
    Laplace = D - Adj;
    e       = ones(n,1);
    degree  = sum(sum(Adj,2));
    L_p     = (Laplace-((e*e')./n))^-1+((e*e')./n);
    ACT     = zeros(n);
    for i = 1 : n
        AI = L_p(i,:);
        for j = 1 : n
            BI = L_p(j,:);
            CI = AI - BI;
            ACT(i,j) = degree*(CI(:,i)-CI(:,j));
        end
    end
end

function C_Community = Decoding(PopDec,ACT,M)
    %% Decoding
    PopDec      = logical(PopDec);
    Community   = 1 : M;
    IndexMatrix = 1 : M;
    if sum(PopDec) > 1
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
end

function [U_new,IndexMatrix] = Myfunction_random(Community,PopDec,IndexMatrix,ACT)
    [~,n] = size(Community);
    tmp   = ACT(PopDec,~PopDec);
    U_new = tmp./(ones(n,1)*sum(tmp));
    IndexMatrix = IndexMatrix(~PopDec);
end

function cf = KKM(Adj,C)
% Calculate the kernel k-means

    cf = 0;
    ec = 0;
    numVar  = size(Adj,1);
    clu_num = length(C);
    for i = 1 : clu_num
        s_index = C{i};
        s = Adj(s_index,s_index);
        s_cardinality = length(s_index);
        if s_cardinality > 0
            kins_sum = 0;
            for j = 1 : s_cardinality
                kins = sum(s(j,:));
                kins_sum = kins_sum + kins;
                ec = kins_sum;
            end
            cf_s = ec*1.0/(s_cardinality);
            cf   = cf + cf_s;
        end
    end
    cf = 2*(numVar-length(C)) - cf;
end

function cs = RC(Adj,C)
% Calculate the ratio cut

    cs = 0;
    de = 0;
    clu_num = length(C);
    for i = 1 : clu_num
        s_index = C{i};
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

function Q = Modularity(Adj,C)
% Calculate the modularity

    Q = 0;
    M = sum(sum(Adj))/2;
    for i = 1 : length(C)
        Q = Q + sum(sum(Adj(C{i},C{i})))/2/M - (sum(sum(Adj(C{i},:)))/2/M)^2;
    end
end
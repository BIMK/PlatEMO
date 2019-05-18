classdef Sparse_CN < PROBLEM
% <problem> <Sparse MOP>
% The critical node detection problem
% dataNo --- 1 --- Number of dataset

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

% The datasets are taken from the Pajek datasets in
% http://vlado.fmf.uni-lj.si/pub/networks/data/default.htm
% No.   Name            Nodes   Edges
% 1     Movies           102     192
% 2     GD99             234     154
% 3     GD01             311     640
% 4     GD97             452     460

    properties(Access = private)
        G;  % Adjacency matrix of the graph
    end
    methods
        %% Initialization
        function obj = Sparse_CN()
            % Load data
            dataNo = obj.Global.ParameterSet(1);
            str    = {'Movies','GD99','GD01','GD97'};
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_CN.mat'),'Dataset');
            obj.G = Dataset.(str{dataNo});
            % Parameter setting
            obj.Global.M        = 2;
            obj.Global.D        = length(obj.G);
            obj.Global.encoding = 'binary';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec = logical(PopDec);
            PopObj = zeros(size(PopDec,1),2);
            for i = 1 : size(PopDec,1)
                PopObj(i,1) = mean(PopDec(i,:));
                PopObj(i,2) = PairwiseConnectivity(obj.G(~PopDec(i,:),~PopDec(i,:)));
            end
            PopObj(:,2) = PopObj(:,2)/PairwiseConnectivity(obj.G);
        end
    end
end

function f = PairwiseConnectivity(G)
    f = 0;
    remain = true(1,length(G));
    while any(remain)
        c = false(1,length(G));
        c(find(remain,1)) = true;
        while true
            c1 = any(G(c,:),1) & remain;
            if all(c==c1)
                break;
            else
                c = c1;
            end
        end
        remain = remain & ~c;
        f = f + sum(c)*(sum(c)-1)/2;
    end
end
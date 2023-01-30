classdef Sparse_CN < PROBLEM
% <multi> <binary> <large/none> <expensive/none> <sparse/none>
% The critical node detection problem
% dataNo --- 1 --- Number of dataset

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zhang, C. Wang, and Y. Jin, An evolutionary algorithm for
% large-scale sparse multi-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2020, 24(2): 380-393.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
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
        A;  % Adjacency matrix of the graph
        G;  % The graph object
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Load data
            dataNo = obj.ParameterSet(1);
            str    = {'Movies','GD99','GD01','GD97'};
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_CN.mat'),'Dataset');
            obj.A = Dataset.(str{dataNo});
            obj.G = graph(obj.A&~logical(eye(length(obj.A))));
            % Parameter setting
            obj.M        = 2;
            obj.D        = length(obj.A);
            obj.encoding = 4 + zeros(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec = logical(PopDec);
            PopObj = zeros(size(PopDec,1),2);
            for i = 1 : size(PopDec,1)
                PopObj(i,1) = mean(PopDec(i,:));
                PopObj(i,2) = PairwiseConnectivity(obj.A(~PopDec(i,:),~PopDec(i,:)));
            end
            PopObj(:,2) = PopObj(:,2)/PairwiseConnectivity(obj.A);
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            [~,best] = min(sum(Population.objs,2));
            h = plot(Draw([]),obj.G,'MarkerSize',6,'EdgeColor',[.3 .3 .3]);
            highlight(h,logical(Population(best).dec),'NodeColor','k');
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            Draw(Population.objs,{'Ratio of deleted nodes','Ratio of pairwise connectivity',[]});
        end
    end
end

function f = PairwiseConnectivity(A)
    f = 0;
    remain = true(1,length(A));
    while any(remain)
        c = false(1,length(A));
        c(find(remain,1)) = true;
        while true
            c1 = any(A(c,:),1) & remain;
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
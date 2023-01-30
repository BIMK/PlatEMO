classdef CommunityDetection < PROBLEM
% <single> <label> <large/none> <expensive/none>
% The community detection problem with label based encoding
% dataNo --- 1 --- Number of dataset

%------------------------------- Reference --------------------------------
% Y. Tian, S. Yang, and X. Zhang, An evolutionary multiobjective
% optimization based fuzzy method for overlapping community detection, IEEE
% Transactions on Fuzzy Systems, 2020, 28(11): 2841-2855.
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
        G;      % The graph object
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
            obj.G   = graph(obj.Adj);
            % Parameter setting
            obj.M = 1;
            obj.D = size(obj.Adj,2);
            obj.lower    = 1     + zeros(1,obj.D);
            obj.upper    = obj.D + zeros(1,obj.D);
            obj.encoding = 3     + zeros(1,obj.D);
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj,PopDec)
            for i = 1 : size(PopDec,1)
                P = zeros(1,obj.D);
                while ~all(P)
                    x = find(~P,1);
                    P(PopDec(i,:)==PopDec(i,x)) = max(P) + 1;
                end
                PopDec(i,:) = P;
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopObj,1)
                PopObj(i) = 1 - Modularity(obj.Adj,Decoding(PopDec(i,:)));
            end
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            [~,best] = min(Population.objs);
            C = Decoding(Population(best).dec);
            h = plot(Draw([]),obj.G,'MarkerSize',6,'EdgeColor',[.3 .3 .3]);
            tempStream = RandStream('mlfg6331_64','Seed',2);
            for i = 1 : length(C)
                highlight(h,C{i},'NodeColor',rand(tempStream,1,3));
            end
        end
    end
end

function Community = Decoding(Dec)
    Community = {};
    while any(Dec)
        current      = find(Dec==Dec(find(Dec,1)));
        Community    = [Community,current];
        Dec(current) = 0;
    end
end

function Q = Modularity(Adj,C)
    Q = 0;
    M = sum(sum(Adj))/2;
    for i = 1 : length(C)
        Q = Q + sum(sum(Adj(C{i},C{i})))/2/M - (sum(sum(Adj(C{i},:)))/2/M)^2;
    end
end
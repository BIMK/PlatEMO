classdef MOKP < PROBLEM
% <problem> <Combinatorial MOP>
% The multi-objective knapsack problem

%------------------------------- Reference --------------------------------
% E. Zitzler and L. Thiele, Multiobjective evolutionary algorithms: A
% comparative case study and the strength Pareto approach, IEEE
% Transactions on Evolutionary Computation, 1999, 3(4): 257-271.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        P;	% Profit of each item according to each knapsack
        W;  % Weight of each item according to each knapsack
    end
    methods
        %% Initialization
        function obj = MOKP()
            % Parameter setting
            if isempty(obj.Global.M)
                obj.Global.M = 2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 250;
            end
            obj.Global.encoding = 'binary';
            % Randomly generate profits and weights
            file = sprintf('MOKP-M%d-D%d.mat',obj.Global.M,obj.Global.D);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'P','W');
            else
                P = randi([10,100],obj.Global.M,obj.Global.D);
                W = randi([10,100],obj.Global.M,obj.Global.D);
                save(file,'P','W');
            end
            obj.P = P;
            obj.W = W;
        end
        %% Repair infeasible solutions
        function PopDec = CalDec(obj,PopDec)
            C = sum(obj.W,2)/2;
            [~,rank] = sort(max(obj.P./obj.W));
            for i = 1 : size(PopDec,1)
                while any(obj.W*PopDec(i,:)'>C)
                    k = find(PopDec(i,rank),1);
                    PopDec(i,rank(k)) = 0;
                end
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = repmat(sum(obj.P,2)',size(PopDec,1),1) - PopDec*obj.P';
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
            P = sum(obj.P,2)';
        end
        %% Draw special figure
        function Draw(obj,PopDec)
            cla; Draw(PopDec*obj.P');
        end
    end
end
classdef PROBLEM < handle
%PROBLEM - The superclass of all the problems.
%
%   This is the superclass of all the problems. This class cannot be
%   instantiated.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        Global; % The current GLOBAL object
    end
    methods(Access = protected)
        %% Constructor
        function obj = PROBLEM()
            obj.Global = GLOBAL.GetObj();
        end
    end
    methods
        %% Generate initial population
        function PopDec = Init(obj,N)
            switch obj.Global.encoding
                case 'binary'
                    PopDec = randi([0,1],N,obj.Global.D);
                case 'permutation'
                    [~,PopDec] = sort(rand(N,obj.Global.D),2);
                otherwise
                    PopDec = unifrnd(repmat(obj.Global.lower,N,1),repmat(obj.Global.upper,N,1));
            end
        end
        %% Repair infeasible solutions
        function PopDec = CalDec(obj,PopDec)
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj(:,1) = PopDec(:,1)   + sum(PopDec(:,2:end),2);
            PopObj(:,2) = 1-PopDec(:,1) + sum(PopDec(:,2:end),2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon = zeros(size(PopDec,1),1);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = ones(1,obj.Global.M);
        end
        %% Draw special figure
        function Draw(obj,PopDec)
        end
    end
end
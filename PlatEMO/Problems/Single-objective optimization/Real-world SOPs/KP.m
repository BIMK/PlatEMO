classdef KP < PROBLEM
% <single> <binary> <large/none> <constrained>
% The knapsack problem

%------------------------------- Reference --------------------------------
% E. Zitzler and L. Thiele, Multiobjective evolutionary algorithms: A
% comparative case study and the strength Pareto approach, IEEE
% Transactions on Evolutionary Computation, 1999, 3(4): 257-271.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        P;	% Profit of each item
        W;  % Weight of each item
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            obj.M = 1;
            if isempty(obj.D); obj.D = 250; end
            obj.encoding = 4 + zeros(1,obj.D);
            % Randomly generate profits and weights
            file = sprintf('KP-D%d.mat',obj.D);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'P','W');
            else
                P = randi([10,100],1,obj.D);
                W = randi([10,100],1,obj.D);
                save(file,'P','W');
            end
            obj.P = P;
            obj.W = W;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = repmat(sum(obj.P,2)',size(PopDec,1),1) - PopDec*obj.P';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon = PopDec*obj.W' - sum(obj.W)/2;
        end
    end
end
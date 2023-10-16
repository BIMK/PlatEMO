classdef LIRCMOP1 < PROBLEM
% <multi> <real> <large/none> <constrained>
% Constrained benchmark MOP with large infeasible regions

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Huang, Y. Fang, Y. You, J. Mo, C. Wei, and E.
% Goodman, An improved epsilon constraint-handling method in MOEA/D for
% CMOPs with large infeasible regions, Soft Computing, 2019, 23:
% 12491-12510.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenji Li
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D); obj.D = 30; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            x_odd       = X(:,3:2:end);
            x_even      = X(:,2:2:end);
            g_1         = sum((x_odd - sin(0.5 * pi * X(:,1))).^2,2);
            g_2         = sum((x_even - cos(0.5 * pi * X(:,1))).^2,2);
            PopObj(:,1) = X(:,1) + g_1;
            PopObj(:,2) = 1 - X(:,1) .^ 2 + g_2;
            PopCon(:,1) = (0.5 - g_1).*(0.51 - g_1);
            PopCon(:,2) = (0.5 - g_2).*(0.51 - g_2);
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - R(:,1).^2;
            R      = R + 0.5;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end
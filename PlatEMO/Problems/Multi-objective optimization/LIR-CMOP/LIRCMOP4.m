classdef LIRCMOP4 < PROBLEM
% <multi> <real> <large/none> <constrained>
% Constrained benchmark MOP with large infeasible regions

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Huang, Y. Fang, Y. You, J. Mo, C. Wei, and E.
% Goodman, An improved epsilon constraint-handling method in MOEA/D for
% CMOPs with large infeasible regions, Soft Computing, 2019, 23:
% 12491-12510.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
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
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            x_odd       = X(:,3:2:end);
            x_even      = X(:,2:2:end);
            len_odd     = size(x_odd,2); 
            len_even    = size(x_even,2);
            g_1         = sum((x_odd - repmat(X(:,1),1,len_odd)).^2,2);
            g_2         = sum((x_even - repmat(X(:,1),1,len_even)).^2,2);  
            PopObj(:,1) = X(:,1) + g_1;
            PopObj(:,2) = 1 - sqrt(X(:,1)) + g_2;
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            x_odd       = X(:,3:2:end);
            x_even      = X(:,2:2:end);
            len_odd     = size(x_odd,2); 
            len_even    = size(x_even,2);
            g_1         = sum((x_odd - repmat(X(:,1),1,len_odd)).^2,2);
            g_2         = sum((x_even - repmat(X(:,1),1,len_even)).^2,2);  
            PopCon(:,1) = (0.5 - g_1).*(0.51 - g_1);
            PopCon(:,2) = (0.5 - g_2).*(0.51 - g_2);
            PopCon(:,3) = 0.5 - sin(20 * pi * X(:,1));
        end 
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - sqrt(R(:,1));
            R(sin(20*pi*R(:,1))<0.5,:) = [];
            R      = R + 0.5;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R(:,1) = linspace(0,1,100)';
            R(:,2) = 1 - sqrt(R(:,1));
            R(sin(20*pi*R(:,1))<0.5,:) = nan;
            R      = R + 0.5;
        end
    end
end
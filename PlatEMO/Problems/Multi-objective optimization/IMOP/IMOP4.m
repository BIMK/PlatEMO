classdef IMOP4 < PROBLEM
% <multi> <real> <expensive/none>
% Benchmark MOP with irregular Pareto front
% a1 --- 0.05 --- Parameter a1
% K  ---    5 --- Parameter K

%------------------------------- Reference --------------------------------
% Y. Tian, R. Cheng, X. Zhang, M. Li, and Y. Jin, Diversity assessment of
% multi-objective evolutionary algorithms: Performance metric and benchmark
% problems, IEEE Computational Intelligence Magazine, 2019, 14(3): 61-74.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        a1 = 0.05;  % Parameter a1
        K  = 5;     % Parameter K
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.a1,obj.K] = obj.ParameterSet(0.05,5);
            obj.M = 3;
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            y1 = mean(PopDec(:,1:obj.K),2).^obj.a1;
            g  = sum((PopDec(:,obj.K+1:end)-0.5).^2,2);
            PopObj(:,1) = (1+g).*(y1);
            PopObj(:,2) = (1+g).*(y1+sin(10*pi*y1)/10);
            PopObj(:,3) = (1+g).*(1-y1);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = R(:,1) + sin(10*pi*R(:,1))/10;
            R(:,3) = 1 - R(:,1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end
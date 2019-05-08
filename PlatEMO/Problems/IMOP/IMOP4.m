classdef IMOP4 < PROBLEM
% <problem> <IMOP>
% Benchmark MOP with irregular Pareto front
% a1 --- 0.05 --- Parameter a1
% K  ---    5 --- Parameter K

%------------------------------- Reference --------------------------------
% Y. Tian, R. Cheng, X. Zhang, M. Li, and Y. Jin, Diversity assessment of
% multi-objective evolutionary algorithms: Performance metric and benchmark
% problems, IEEE Computational Intelligence Magazine, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
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
        %% Initialization
        function obj = IMOP4()
            [obj.a1,obj.K] = obj.Global.ParameterSet(0.05,5);
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            y1 = mean(PopDec(:,1:obj.K),2).^obj.a1;
            g  = sum((PopDec(:,obj.K+1:end)-0.5).^2,2);
            PopObj(:,1) = (1+g).*(y1);
            PopObj(:,2) = (1+g).*(y1+sin(10*pi*y1)/10);
            PopObj(:,3) = (1+g).*(1-y1);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = 0 : 1/(N-1) : 1;
            P(:,2) = P(:,1) + sin(10*pi*P(:,1))/10;
            P(:,3) = 1 - P(:,1);
        end
    end
end
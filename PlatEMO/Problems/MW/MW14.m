classdef MW14 < PROBLEM
% <problem> <MW>
% Constrained benchmark MOP proposed by Ma and Wang

%------------------------------- Reference --------------------------------
% Z. Ma and Y. Wang, Evolutionary constrained multiobjective optimization:
% Test suite construction and performance comparisons. IEEE Transactions on
% Evolutionary Computation, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function obj = MW14()
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 15;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            X = 1.5*X;
            M = obj.Global.M;
            g = sum(2*(X(:,M:end) + (X(:,M-1:end-1) - 0.5).^2 - 1).^2,2);
            PopObj(:,1:M-1) = X(:,1:M-1);
            PopObj(:,M)     = ((1+g)/(M-1)).*sum(6 - exp(PopObj(:,1:M-1)) - 1.5*sin(1.1*pi*PopObj(:,1:M-1).^2),2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            a      = 1 + PopObj(:,1:obj.Global.M-1) + 0.5*PopObj(:,1:obj.Global.M-1).^2 + 1.5*sin(1.1*pi*PopObj(:,1:obj.Global.M-1).^2);
            PopCon = PopObj(:,obj.Global.M) - 1/(obj.Global.M-1)*sum(6.1 - a,2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            interval     = [0,0.731000,1.331000,1.500000];
            median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
            X            = ReplicatePoint(N,obj.Global.M-1);
            X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
            X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
            P            = [X,1/(obj.Global.M-1)*sum(6 - exp(X) - 1.5*sin(1.1*pi*X.^2),2)];
        end
    end
end

function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end
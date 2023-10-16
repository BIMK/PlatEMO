classdef MW14 < PROBLEM
% <multi/many> <real> <large/none> <constrained>
% Constrained benchmark MOP proposed by Ma and Wang

%------------------------------- Reference --------------------------------
% Z. Ma and Y. Wang, Evolutionary constrained multiobjective optimization:
% Test suite construction and performance comparisons. IEEE Transactions on
% Evolutionary Computation, 2019, 23(6): 972-986.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = 15; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X      = varargin{1};
            X      = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            PopDec = X;
            X      = 1.5*X;
            g      = sum(2*(X(:,obj.M:end) + (X(:,obj.M-1:end-1) - 0.5).^2 - 1).^2,2);
            PopObj(:,1:obj.M-1) = X(:,1:obj.M-1);
            PopObj(:,obj.M)     = ((1+g)/(obj.M-1)).*sum(6 - exp(PopObj(:,1:obj.M-1)) - 1.5*sin(1.1*pi*PopObj(:,1:obj.M-1).^2),2);
            a           = 1 + PopObj(:,1:obj.M-1) + 0.5*PopObj(:,1:obj.M-1).^2 + 1.5*sin(1.1*pi*PopObj(:,1:obj.M-1).^2);
            PopCon      = PopObj(:,obj.M) - 1/(obj.M-1)*sum(6.1 - a,2);
            Population  = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            interval     = [0,0.731000,1.331000,1.500000];
            median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
            X            = UniformPoint(N,obj.M-1,'grid');
            X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
            X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
            R            = [X,1/(obj.M-1)*sum(6 - exp(X) - 1.5*sin(1.1*pi*X.^2),2)];
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            if obj.M == 2
                x      = linspace(0,1.5,100)';
                y      = 6 - exp(x) - 1.5*sin(1.1*pi*x.^2);
                nd     = NDSort([x,y],1)==1;
                x(~nd) = nan;
                R      = [x,y];
            elseif obj.M == 3
                [x,y]  = meshgrid(linspace(0,1.5,20));
                z      = 1/2*(12-exp(x)-1.5*sin(1.1*pi*x.^2)-exp(y)-1.5*sin(1.1*pi*y.^2));
                nd     = reshape(NDSort([x(:),y(:),z(:)],1)==1,size(z));
                z(~nd) = nan;
                R      = {x,y,z};
            else
                R = [];
            end
        end
    end
end
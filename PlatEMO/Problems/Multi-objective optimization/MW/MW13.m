classdef MW13 < PROBLEM
% <multi> <real> <large/none> <constrained>
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
            obj.M = 2;
            if isempty(obj.D); obj.D = 15; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            z = 1 - exp(-10*(X(:,obj.M:end) - (repmat(obj.M:obj.D,size(X,1),1) - 1)/obj.D).^2);
            g = 1 + sum((1.5 + (0.1/obj.D)*z.^2 - 1.5*cos(2*pi*z)),2);
            PopObj(:,1) = g.*X(:,1)*1.5;
            PopObj(:,2) = g.*(5 - exp(PopObj(:,1)./g) - abs(0.5*sin(3*pi*PopObj(:,1)./g)));
            PopCon(:,1) = (5 - exp(PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2)).*(5 - (1 + 0.4*PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2));
            PopCon(:,2) = -(5 - (1 + PopObj(:,1) + 0.5*PopObj(:,1).^2) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2)).*(5 - (1 + 0.7*PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2));
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1)  = (0:1.5/(N-1):1.5)';
            R(:,2)  = 5 - exp(R(:,1)) - 0.5*abs(sin(3*pi*R(:,1)));
            c1      = (5-exp(R(:,1))-0.5*sin(3*pi*R(:,1))-R(:,2)).*(5-(1+0.4*R(:,1))-0.5*sin(3*pi*R(:,1))-R(:,2));
            invalid = c1>0;
            while any(invalid)
                R(invalid,:) = R(invalid,:).*1.001;
                c1      = (5-exp(R(:,1))-0.5*sin(3*pi*R(:,1))-R(:,2)).*(5-(1+0.4*R(:,1))-0.5*sin(3*pi*R(:,1))-R(:,2));
                invalid = c1>0;
            end
            R = R(NDSort(R,1)==1,:);
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0,2,400),linspace(0,4.5,400));
            z     = nan(size(x));
            fes1  = (5-exp(x)-0.5*sin(3*pi*x)-y).*(5-(1+0.4*x)-0.5*sin(3*pi*x)-y) <= 0;
            fes2  = -(5-(1+x+0.5*x.^2)-0.5*sin(3*pi*x)-y).*(5-(1+0.7*x)-0.5*sin(3*pi*x)-y) <= 0;
            z(fes1 & fes2 & exp(x)+abs(0.5*sin(3*pi*x))+y>=5) = 0;
            R = {x,y,z};
        end
    end
end
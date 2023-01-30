classdef MW8 < PROBLEM
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
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            z = 1 - exp(-10*(X(:,obj.M:end) - ((obj.M:obj.D) - 1)/obj.D).^2);
            g = sum((1.5 + (0.1/obj.D)*z.^2 - 1.5*cos(2*pi*z)),2);
            PopObj = repmat(1+g,1,obj.M).*flip(cumprod([ones(size(X,1),1),cos(X(:,1:obj.M-1)*pi/2)],2),2).*[ones(size(X,1),1),sin(X(:,obj.M-1:-1:1)*pi/2)];
            l      = asin(PopObj(:,end)./sqrt(sum(PopObj.^2,2)));
            PopCon = sum(PopObj.^2,2) - (1.25 - 0.5*sin(6*l).^2).^2;
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
             R = UniformPoint(N,obj.M);
             R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
             R(1-(1.25 - 0.5*sin(6*asin(R(:,end))).^2).^2>0,:) = [];
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            if obj.M == 2
                [x,y] = meshgrid(linspace(0,1.2,400));
                z     = nan(size(x));
                fes   = x.^2 + y.^2 - (1.25-0.5*sin(6*asin(y./sqrt(x.^2+y.^2))).^2).^2 <= 0;
                z(fes & x.^2+y.^2>=1) = 0;
                R = {x,y,z};
            elseif obj.M == 3
                a = linspace(0,pi/2,40)';
                x = sin(a)*cos(a');
                y = sin(a)*sin(a');
                z = cos(a)*ones(size(a'));
                fes     = 1 - (1.25-0.5*sin(6*asin(z)).^2).^2 <= 0;
                z(~fes) = nan;
                R = {x,y,z};
            else
                R = [];
            end
        end
    end
end
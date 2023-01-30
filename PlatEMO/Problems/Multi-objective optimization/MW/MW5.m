classdef MW5 < PROBLEM
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
            g = 1 + sum(1 - exp(-10*((X(:,obj.M:end).^(obj.D-obj.M)) - 0.5 - (repmat(obj.M:obj.D,size(X,1),1) - 1)/(2*obj.D)).^2),2);
            PopObj(:,1) = g.*X(:,1);
            PopObj(:,2) = g.*sqrt(1 - (PopObj(:,1)./g).^2);
            l1          = atan(PopObj(:,2)./PopObj(:,1));
            l2          = 0.5*pi - 2*abs(l1-0.25*pi);
            PopCon(:,1) = PopObj(:,1).^2 + PopObj(:,2).^2 - (1.7-0.2*sin(2*l1)).^2;
            PopCon(:,2) = (1+0.5*sin(6*l2.^3)).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
            PopCon(:,3) = (1-0.45*sin(6*l2.^3)).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = [0 1;0.3922 0.9199;0.4862 0.8739;0.5490 0.8358;0.5970 0.8023;0.6359 0.7719;0.6686 0.7436;0.6969 0.7174];
            R = [R;flip(R,2)];
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0,1.7,400));
            z     = nan(size(x));
            l1    = atan(y./x);
            l2    = 0.5*pi - 2*abs(l1-0.25*pi);
            fes1  = x.^2 + y.^2 - (1.7-0.2*sin(2*l1)).^2 <= 0;
            fes2  = (1+0.5*sin(6*l2.^3)).^2 - x.^2 - y.^2 <= 0;
            fes3  = (1-0.45*sin(6*l2.^3)).^2 - x.^2 - y.^2 <= 0;
            z(fes1 & fes2 & fes3 & x.^2+y.^2>=1) = 0;
            R = {x,y,z};
        end
    end
end
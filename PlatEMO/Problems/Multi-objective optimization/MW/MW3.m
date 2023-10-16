classdef MW3 < PROBLEM
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
            g = 1 + sum(2*(X(:,obj.M:end) + (X(:,obj.M-1:end-1) - 0.5).^2 - 1).^2,2);
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1 - PopObj(:,1)./g);
            l           = sqrt(2)*PopObj(:,2) - sqrt(2)*PopObj(:,1);
            PopCon(:,1) = sum(PopObj,2) - 1.05 - 0.45*sin(0.75*pi*l).^6;
            PopCon(:,2) = 0.85 - sum(PopObj,2) + 0.3*sin(0.75*pi*l).^2;
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1)  = (0:1/(N-1):1)';
            R(:,2)  = 1 - R(:,1);
            invalid = (0.85-R(:,1)-R(:,2)+0.3*sin(0.75*pi*sqrt(2)*(R(:,2)-R(:,1))).^2) > 0;
            while any(invalid)
                R(invalid,:) = R(invalid,:).*1.001;
                invalid = (0.85-R(:,1)-R(:,2)+0.3*sin(0.75*pi*sqrt(2)*(R(:,2)-R(:,1))).^2) > 0;
            end
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0,1,400),linspace(0,1.5,400));
            z     = nan(size(x));
            fes1  = x + y - 1.05 - 0.45*sin(0.75*pi*(sqrt(2)*y-sqrt(2)*x)).^6 <= 0;
            fes2  = 0.85 - x - y + 0.3*sin(0.75*pi*(sqrt(2)*y-sqrt(2)*x)).^2 <= 0;
            z(fes1 & fes2 & x+y>=1) = 0;
            R = {x,y,z};
        end
    end
end
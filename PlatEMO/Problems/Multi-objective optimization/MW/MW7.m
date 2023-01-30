classdef MW7 < PROBLEM
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
            PopObj(:,1) = g.*X(:,1);
            PopObj(:,2) = g.*sqrt(1 - (PopObj(:,1)./g).^2);
            l           = atan(PopObj(:,2)./PopObj(:,1));
            PopCon(:,1) = PopObj(:,1).^2 + PopObj(:,2).^2 - (1.2+0.4*sin(4*l).^16).^2;
            PopCon(:,2) = (1.15-0.2*sin(4*l).^8).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = (0:1/(N-1):1)';
            R(:,2) = 1 - R(:,1);
            R = R./repmat(sqrt(sum(R.^2,2)),1,2);
            invalid = ((1.15-0.2*sin(4*atan(R(:,2)./R(:,1))).^8).^2-R(:,1).^2-R(:,2).^2) > 0;
            while any(invalid)
                R(invalid,:) = R(invalid,:).*1.001;
                invalid = ((1.15-0.2*sin(4*atan(R(:,2)./R(:,1))).^8).^2-R(:,1).^2-R(:,2).^2) > 0;
            end
            R = R(NDSort(R,1)==1,:);
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0,1.5,400));
            z     = nan(size(x));
            l     = atan(y./x);
            fes1  = x.^2 + y.^2 - (1.2+0.4*sin(4*l).^16).^2 <= 0;
            fes2  = (1.15-0.2*sin(4*l).^8).^2 - x.^2 - y.^2 <= 0;
            z(fes1 & fes2 & x.^2+y.^2>=1) = 0;
            R = {x,y,z};
        end
    end
end
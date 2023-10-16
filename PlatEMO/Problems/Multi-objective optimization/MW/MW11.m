classdef MW11 < PROBLEM
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
            PopObj(:,1) = g.*X(:,1)*sqrt(1.9999);
            PopObj(:,2) = g.*sqrt(2 - (PopObj(:,1)./g).^2);
            PopCon(:,1) = -(3 - PopObj(:,1).^2 - PopObj(:,2)).*(3 - 2*PopObj(:,1).^2 - PopObj(:,2));
            PopCon(:,2) = (3 - 0.625*PopObj(:,1).^2 - PopObj(:,2)).*(3 - 7*PopObj(:,1).^2 - PopObj(:,2));
            PopCon(:,3) = -(1.62 - 0.18*PopObj(:,1).^2 - PopObj(:,2)).*(1.125 - 0.125*PopObj(:,1).^2 - PopObj(:,2));
            PopCon(:,4) = (2.07 - 0.23*PopObj(:,1).^2 - PopObj(:,2)).*(0.63 - 0.07*PopObj(:,1).^2 - PopObj(:,2));
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1)  = (0:1/(N-1):1)';
            R(:,2)  = 1 - R(:,1);
            R       = R./repmat(sqrt(sum(R.^2,2)/2),1,2);
            c1      = (3 - R(:,1).^2 - R(:,2)).*(3 - 2*R(:,1).^2 - R(:,2));
            c2      = (3 - 0.625*R(:,1).^2 - R(:,2)).*(3 - 7*R(:,1).^2 - R(:,2));
            c3      = (1.62 - 0.18*R(:,1).^2 - R(:,2)).*(1.125 - 0.125*R(:,1).^2 - R(:,2));
            c4      = (2.07 - 0.23*R(:,1).^2 - R(:,2)).*(0.63 - 0.07*R(:,1).^2 - R(:,2));
            invalid = c1<0 | c2>0 | c3<0 | c4>0;
            while any(invalid)
                R(invalid,:) = R(invalid,:).*1.001;
                R(any(R>2.2,2),:) = [];
                c1      = (3 - R(:,1).^2 - R(:,2)).*(3 - 2*R(:,1).^2 - R(:,2));
                c2      = (3 - 0.625*R(:,1).^2 - R(:,2)).*(3 - 7*R(:,1).^2 - R(:,2));
                c3      = (1.62 - 0.18*R(:,1).^2 - R(:,2)).*(1.125 - 0.125*R(:,1).^2 - R(:,2));
                c4      = (2.07 - 0.23*R(:,1).^2 - R(:,2)).*(0.63 - 0.07*R(:,1).^2 - R(:,2));
                invalid = c1<0 | c2>0 | c3<0 | c4>0;
            end
            R = [R;1,1];
            R = R(NDSort(R,1)==1,:);
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0,2.1,400));
            z     = nan(size(x));
            fes1  = -(3-x.^2-y).*(3-2*x.^2-y) <= 0;
            fes2  = (3-0.625*x.^2-y).*(3-7*x.^2-y) <= 0;
            fes3  = -(1.62-0.18*x.^2-y).*(1.125-0.125*x.^2-y) <= 0;
            fes4  = (2.07-0.23*x.^2-y).*(0.63-0.07*x.^2-y) <= 0;
            z(fes1 & fes2 & fes3 & fes4 & x.^2+y.^2>=2) = 0;
            R = {x,y,z};
        end
    end
end
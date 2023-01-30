classdef MW12 < PROBLEM
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
            PopObj(:,2) = g.*(0.85 - 0.8*(PopObj(:,1)./g) - 0.08*abs(sin(3.2*pi*(PopObj(:,1)./g))));
            PopCon(:,1) = (1 - 0.8*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2) - PopObj(:,1)/1.5))).*(1.8 - 1.125*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2)/1.8 - PopObj(:,1)/1.6)));
            PopCon(:,2) = -(1 - 0.625*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2) - PopObj(:,1)/1.6))).*(1.4 - 0.875*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2)/1.4 - PopObj(:,1)/1.6)));
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1)  = (0:1/(N-1):1)';
            R(:,2)  = 0.85 - 0.8*R(:,1) - 0.08*abs(sin(3.2*pi*R(:,1)));
            c1      = (1-0.8*R(:,1)-R(:,2)+0.08*sin(2*pi*(R(:,2)-R(:,1)/1.5))).*(1.8-1.125*R(:,1)-R(:,2)+0.08*sin(2*pi*(R(:,2)/1.8-R(:,1)/1.6)));
            invalid = c1>0;
            while any(invalid)
                R(invalid,:) = R(invalid,:).*1.001;
                c1      = (1-0.8*R(:,1)-R(:,2)+0.08*sin(2*pi*(R(:,2)-R(:,1)/1.5))).*(1.8-1.125*R(:,1)-R(:,2)+0.08*sin(2*pi*(R(:,2)/1.8-R(:,1)/1.6)));
                invalid = c1>0;
            end
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0,2,400));
            z     = nan(size(x));
            fes1  = (1-0.8*x-y+0.08*sin(2*pi*(y-x/1.5))).*(1.8-1.125*x-y+0.08*sin(2*pi*(y/1.8-x/1.6))) <= 0;
            fes2  = -(1-0.625*x-y+0.08*sin(2*pi*(y-x/1.6))).*(1.4-0.875*x-y+0.08*sin(2*pi*(y/1.4-x/1.6))) <= 0;
            z(fes1 & fes2 & 0.8*x+0.08*abs(sin(3.2*pi*x))+y>=0.85) = 0;
            R = {x,y,z};
        end
    end
end
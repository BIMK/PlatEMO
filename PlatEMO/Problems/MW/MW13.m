classdef MW13 < PROBLEM
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
        function obj = MW13()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 15;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,D] = size(X);
            M     = obj.Global.M;
            z     = 1 - exp(-10*(X(:,M:end) - (repmat(M:D,N,1) - 1)/D).^2);
            g     = 1 + sum((1.5 + (0.1/D)*z.^2 - 1.5*cos(2*pi*z)),2);
            PopObj(:,1) = g.*X(:,1)*1.5;
            PopObj(:,2) = g.*(5 - exp(PopObj(:,1)./g) - abs(0.5*sin(3*pi*PopObj(:,1)./g)));
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            PopCon(:,1) = (5 - exp(PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2)).*(5 - (1 + 0.4*PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2));
            PopCon(:,2) = -(5 - (1 + PopObj(:,1) + 0.5*PopObj(:,1).^2) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2)).*(5 - (1 + 0.7*PopObj(:,1)) - 0.5*sin(3*pi*PopObj(:,1)) - PopObj(:,2));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1.5/(N-1):1.5)';
            P(:,2) = 5 - exp(P(:,1)) - 0.5*abs(sin(3*pi*P(:,1)));
            c1      = (5-exp(P(:,1))-0.5*sin(3*pi*P(:,1))-P(:,2)).*(5-(1+0.4*P(:,1))-0.5*sin(3*pi*P(:,1))-P(:,2));
            invalid = c1>0;
            while any(invalid)
                P(invalid,:) = P(invalid,:).*1.001;
                c1      = (5-exp(P(:,1))-0.5*sin(3*pi*P(:,1))-P(:,2)).*(5-(1+0.4*P(:,1))-0.5*sin(3*pi*P(:,1))-P(:,2));
                invalid = c1>0;
            end
            P = P(NDSort(P,1)==1,:);
        end
    end
end
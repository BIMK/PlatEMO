classdef MW3 < PROBLEM
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
        function obj = MW3()
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
            M      = obj.Global.M;
            g      = 1 + sum(2*(X(:,M:end) + (X(:,M-1:end-1) - 0.5).^2 - 1).^2,2);
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1 - PopObj(:,1)./g);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            l      = sqrt(2)*PopObj(:,2) - sqrt(2)*PopObj(:,1);
            PopCon(:,1) = sum(PopObj,2) - 1.05 - 0.45*sin(0.75*pi*l).^6;
            PopCon(:,2) = 0.85 - sum(PopObj,2) + 0.3*sin(0.75*pi*l).^2;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1)  = (0:1/(N-1):1)';
            P(:,2)  = 1 - P(:,1);
            invalid = (0.85-P(:,1)-P(:,2)+0.3*sin(0.75*pi*sqrt(2)*(P(:,2)-P(:,1))).^2) > 0;
            while any(invalid)
                P(invalid,:) = P(invalid,:).*1.001;
                invalid = (0.85-P(:,1)-P(:,2)+0.3*sin(0.75*pi*sqrt(2)*(P(:,2)-P(:,1))).^2) > 0;
            end
        end
    end
end
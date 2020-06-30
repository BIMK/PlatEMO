classdef MW7 < PROBLEM
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
        function obj = MW7()
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
            M = obj.Global.M;
            g = 1 + sum(2*(X(:,M:end) + (X(:,M-1:end-1) - 0.5).^2 - 1).^2,2);
            PopObj(:,1) = g.*X(:,1);
            PopObj(:,2) = g.*sqrt(1 - (PopObj(:,1)./g).^2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            l      = atan(PopObj(:,2)./PopObj(:,1));
            PopCon(:,1) = PopObj(:,1).^2 + PopObj(:,2).^2 - (1.2+0.4*sin(4*l).^16).^2;
            PopCon(:,2) = (1.15-0.2*sin(4*l).^8).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1);
            P = P./repmat(sqrt(sum(P.^2,2)),1,2);
            invalid = ((1.15-0.2*sin(4*atan(P(:,2)./P(:,1))).^8).^2-P(:,1).^2-P(:,2).^2) > 0;
            while any(invalid)
                P(invalid,:) = P(invalid,:).*1.001;
                invalid = ((1.15-0.2*sin(4*atan(P(:,2)./P(:,1))).^8).^2-P(:,1).^2-P(:,2).^2) > 0;
            end
            P = P(NDSort(P,1)==1,:);
        end
    end
end
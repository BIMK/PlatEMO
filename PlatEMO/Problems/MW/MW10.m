classdef MW10 < PROBLEM
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
        function obj = MW10()
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
            PopObj(:,1) = g.*(X(:,1).^D);
            PopObj(:,2) = g.*(1 - (PopObj(:,1)./g).^2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            PopCon(:,1) = -(2 - 4*PopObj(:,1).^2 - PopObj(:,2)).*(2 - 8*PopObj(:,1).^2 - PopObj(:,2));
            PopCon(:,2) = (2 - 2*PopObj(:,1).^2 - PopObj(:,2)).*(2 - 16*PopObj(:,1).^2 - PopObj(:,2));
            PopCon(:,3) = (1 - PopObj(:,1).^2 - PopObj(:,2)).*(1.2 - 1.2*PopObj(:,1).^2 - PopObj(:,2));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1).^2;
            c1      = (2 - 4*P(:,1).^2 - P(:,2)).*(2 - 8*P(:,1).^2 - P(:,2));
            c2      = (2 - 2*P(:,1).^2 - P(:,2)).*(2 - 16*P(:,1).^2 - P(:,2));
            c3      = (1 - P(:,1).^2 - P(:,2)).*(1.2 - 1.2*P(:,1).^2 - P(:,2));
            invalid = c1<0 | c2>0 | c3>0;
            while any(invalid)
                P(invalid,:) = P(invalid,:).*1.001;
                P(any(P>1.3,2),:) = [];
                c1      = (2 - 4*P(:,1).^2 - P(:,2)).*(2 - 8*P(:,1).^2 - P(:,2));
                c2      = (2 - 2*P(:,1).^2 - P(:,2)).*(2 - 16*P(:,1).^2 - P(:,2));
                c3      = (1 - P(:,1).^2 - P(:,2)).*(1.2 - 1.2*P(:,1).^2 - P(:,2));
                invalid = c1<0 | c2>0 | c3>0;
            end
            P = P(NDSort(P,1)==1,:);
        end
    end
end
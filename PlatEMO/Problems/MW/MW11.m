classdef MW11 < PROBLEM
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
        function obj = MW11()
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
            PopObj(:,1) = g.*X(:,1)*sqrt(1.9999);
            PopObj(:,2) = g.*sqrt(2 - (PopObj(:,1)./g).^2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            PopCon(:,1) = -(3 - PopObj(:,1).^2 - PopObj(:,2)).*(3 - 2*PopObj(:,1).^2 - PopObj(:,2));
            PopCon(:,2) = (3 - 0.625*PopObj(:,1).^2 - PopObj(:,2)).*(3 - 7*PopObj(:,1).^2 - PopObj(:,2));
            PopCon(:,3) = -(1.62 - 0.18*PopObj(:,1).^2 - PopObj(:,2)).*(1.125 - 0.125*PopObj(:,1).^2 - PopObj(:,2));
            PopCon(:,4) = (2.07 - 0.23*PopObj(:,1).^2 - PopObj(:,2)).*(0.63 - 0.07*PopObj(:,1).^2 - PopObj(:,2));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1);
            P = P./repmat(sqrt(sum(P.^2,2)/2),1,2);
            c1      = (3 - P(:,1).^2 - P(:,2)).*(3 - 2*P(:,1).^2 - P(:,2));
            c2      = (3 - 0.625*P(:,1).^2 - P(:,2)).*(3 - 7*P(:,1).^2 - P(:,2));
            c3      = (1.62 - 0.18*P(:,1).^2 - P(:,2)).*(1.125 - 0.125*P(:,1).^2 - P(:,2));
            c4      = (2.07 - 0.23*P(:,1).^2 - P(:,2)).*(0.63 - 0.07*P(:,1).^2 - P(:,2));
            invalid = c1<0 | c2>0 | c3<0 | c4>0;
            while any(invalid)
                P(invalid,:) = P(invalid,:).*1.001;
                P(any(P>2.2,2),:) = [];
                c1      = (3 - P(:,1).^2 - P(:,2)).*(3 - 2*P(:,1).^2 - P(:,2));
                c2      = (3 - 0.625*P(:,1).^2 - P(:,2)).*(3 - 7*P(:,1).^2 - P(:,2));
                c3      = (1.62 - 0.18*P(:,1).^2 - P(:,2)).*(1.125 - 0.125*P(:,1).^2 - P(:,2));
                c4      = (2.07 - 0.23*P(:,1).^2 - P(:,2)).*(0.63 - 0.07*P(:,1).^2 - P(:,2));
                invalid = c1<0 | c2>0 | c3<0 | c4>0;
            end
            P = [P;1,1];
            P = P(NDSort(P,1)==1,:);
        end
    end
end
classdef MW6 < PROBLEM
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
        function obj = MW6()
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
            PopObj(:,1) = g.*X(:,1)*1.0999;
            PopObj(:,2) = g.*sqrt(1.1*1.1 - (PopObj(:,1)./g).^2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            l      = cos(6*atan(PopObj(:,2)./PopObj(:,1)).^4).^10;
            PopCon = (PopObj(:,1)./(1+0.15*l)).^2 + (PopObj(:,2)./(1+0.75*l)).^2 - 1;
       end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1);
            P = P./repmat(sqrt(sum(P.^2,2)/1.21),1,2);
            l = cos(6*atan(P(:,2)./P(:,1)).^4).^10;
            c = 1 - (P(:,1)./(1+0.15*l)).^2 - (P(:,2)./(1+0.75*l)).^2;
            P(c<0,:) = [];
        end
    end
end

function answer = LA(A,B,C,D,thera)
    t = thera.^C;
    answer = A*cos(B*t).^D; 
end
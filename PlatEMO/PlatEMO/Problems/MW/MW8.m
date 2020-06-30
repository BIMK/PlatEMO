classdef MW8 < PROBLEM
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
        function obj = MW8()
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 15;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
             [N,D]  = size(X);
             M      = obj.Global.M;
             z      = 1 - exp(-10*(X(:,M:end) - ((M:D) - 1)/D).^2);
             g      = sum((1.5 + (0.1/D)*z.^2 - 1.5*cos(2*pi*z)),2);
             PopObj = repmat(1+g,1,M).*flip(cumprod([ones(N,1),cos(X(:,1:M-1)*pi/2)],2),2).*[ones(N,1),sin(X(:,M-1:-1:1)*pi/2)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            l      = asin(PopObj(:,end)./sqrt(sum(PopObj.^2,2)));
            PopCon = sum(PopObj.^2,2) - (1.25 - 0.5*sin(6*l).^2).^2;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
             P = UniformPoint(N,obj.Global.M);
             P = P./repmat(sqrt(sum(P.^2,2)),1,obj.Global.M);
             P(-((1.25 - 0.5*sin(6*asin(P(:,end))).^2).^2 - sum(P.^2,2))>0,:) = [];
        end
    end
end
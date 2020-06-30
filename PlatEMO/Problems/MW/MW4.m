classdef MW4 < PROBLEM
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
        function obj = MW4()
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
             g      = sum(1 - exp(-10*((X(:,M:end).^(D-M)) - 0.5 - ((M:D) - 1)/(2*D)).^2),2);
             PopObj = repmat(1+g,1,M).*flip(cumprod([ones(N,1),X(:,1:M-1)],2),2).*[ones(N,1),1-X(:,M-1:-1:1)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            l      = PopObj(:,end) - sum(PopObj(:,1:(end-1)),2);
            PopCon = sum(PopObj,2) - (1 + 0.4*sin(2.5*pi*l).^8);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,obj.Global.M);
            l = P(:,end) - sum(P(:,1:end-1),2);
            c = (1+0.4*sin(2.5*pi*l).^8) - sum(P,2);
            P(c<0,:) = [];
        end
    end
end
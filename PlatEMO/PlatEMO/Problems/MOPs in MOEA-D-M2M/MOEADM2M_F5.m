classdef MOEADM2M_F5 < PROBLEM
% <problem> <MOEADM2M>
% Benchmark MOP for MOEA/D-M2M

%------------------------------- Reference --------------------------------
% H. Liu, F. Gu, and Q. Zhang, Decomposition of a multiobjective
% optimization problem into a number of simple multiobjective subproblems,
% IEEE Transactions on Evolutionary Computation, 2014, 18(3): 450-455.
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
        function obj = MOEADM2M_F5()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            t = X(:,2:end) - repmat(sin(pi/2*X(:,1)),1,size(X,2)-1);
            g = 2*abs(cos(pi*X(:,1))).*sum(-0.9*t.^2+abs(t).^0.6,2);
            PopObj(:,1) = (1+g).*X(:,1);
            PopObj(:,2) = (1+g).*(1-sqrt(X(:,1)));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - sqrt(P(:,1));
        end
    end
end
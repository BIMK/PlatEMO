classdef DOC8 < PROBLEM
% <multi> <real> <constrained>
% Benchmark MOP with constraints in both decision and objective spaces

%------------------------------- Reference --------------------------------
% Z. Liu and Y. Wang, Handling constrained multiobjective optimization
% problems with constraints in both the decision and objective spaces. IEEE
% Transactions on Evolutionary Computation, 2019, 23(5): 870-884.
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
            obj.M = 3;
            obj.D = 10;
            obj.lower    = [0 0 500 1000 5000 100 100 100 100 100];
            obj.upper    = [1 1 1000 2000 6000 500 500 500 500 500];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            g_temp = X(:, 3) + X(:, 4) + X(:, 5);
            g = g_temp-7049.2480205286 +1;
            PopObj(:,1) = (X(:,1).*X(:,2)).*g;
            PopObj(:,2) = (X(:,1).*(1 - X(:,2))).*g;
            PopObj(:,3) = (1-X(:,1)).*g;
            % Constraints in objective space
            c(:,1) = max( - (PopObj(:,3) - 0.4).*(PopObj(:,3) - 0.6), 0);
            % Constraints in decision space
            c(:,2) = -1 + 0.0025 * (X(:, 6) + X(:, 8));
            c(:,3) = -1 + 0.0025 * (X(:, 7) + X(:, 9) - X(:, 6));
            c(:,4) = -1 + 0.01 * (X(:, 10) - X(:, 7));
            c(:,5) = -X(:, 3).* X(:, 8) + 833.33252 * X(:, 6) + 100 * X(:, 3) - 83333.333;
            c(:,6) = -X(:, 4).* X(:, 9) + 1250 * X(:, 7) + X(:, 4).* X(:, 6) - 1250 * X(:, 6);
            c(:,7) = -X(:, 5).* X(:, 10) + 1250000 + X(:, 5).* X(:, 7) - 2500 * X(:, 7);
            Population  = SOLUTION(X,PopObj,c,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
             R = UniformPoint(N,3);
             R(0.4<R(:,3)&R(:,3)<0.6,:) = [];
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            a = linspace(0,1,20)';
            x = a*a';
            y = a*(1-a');
            z = (1-a)*ones(size(a'));
            z(0.4<z&z<0.6) = nan;
            R = {x,y,z};
        end
    end
end
classdef DOC9 < PROBLEM
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
            obj.D = 11;
            obj.lower    = [0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1];
            obj.upper    = [1 1 10 10 10 10 10 10 10 10 10];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            g_temp = -0.5 * (X(:, 3).* X(:, 6) - X(:, 4).* X(:, 5) + X(:, 5).* X(:, 11) - X(:, 7).* X(:, 11) + X(:, 7).* X(:, 10) - X(:, 8).* X(:, 9));
            g = g_temp + 0.8660254038 +1;
            PopObj(:,1) = cos(0.5*pi*X(:,1)).*cos(0.5*pi*X(:,2)).*g;
            PopObj(:,2) = cos(0.5*pi*X(:,1)).*sin(0.5*pi*X(:,2)).*g;
            PopObj(:,3) = sin(0.5*pi*X(:,1)).*g;
            % Constraints in objective space
            c(:,1) = max( -( PopObj(:,1).^2 + PopObj(:,2).^2 - 1), 0);
            % Constraints in decision space
            c(:,2) = X(:, 5).^2 + X(:, 6).^2 - 1;
            c(:,3) = X(:, 11).^2 - 1;
            c(:,4) = X(:, 7).^2 + X(:, 8).^2 - 1;
            c(:,5) = X(:, 3).^2 + (X(:, 4) - X(:, 11)).^2 - 1;
            c(:,6) = (X(:, 3) - X(:, 7)).^2 + (X(:, 4) - X(:, 8)).^2 - 1;
            c(:,7) = (X(:, 3) - X(:, 9)).^2 + (X(:, 4) - X(:, 10)).^2 - 1;
            c(:,8) = (X(:, 5) - X(:, 7)).^2 + (X(:, 6) - X(:, 8)).^2 - 1;
            c(:,9) = (X(:, 5) - X(:, 9)).^2 + (X(:, 6) - X(:, 10)).^2 - 1;
            c(:,10) = X(:, 9).^2 + (X(:, 10) - X(:, 11)).^2 - 1;
            c(:,11) = X(:, 4).* X(:, 5) - X(:, 3).* X(:, 6);
            c(:,12) = -X(:, 5).* X(:, 11);
            c(:,13) = X(:, 7).* X(:, 11);
            c(:,14) = X(:, 8).* X(:, 9) - X(:, 7).* X(:, 10);
            Population  = SOLUTION(X,PopObj,c,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,2);
            R = R./repmat(sqrt(sum(R.^2,2)),1,2);
            R = [R,zeros(size(R,1),1)];
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end
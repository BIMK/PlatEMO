classdef DOC1 < PROBLEM
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
            obj.M = 2;
            obj.D = 6;
            obj.lower    = [0 78 33 27 27 27];
            obj.upper    = [1 102 45 45 45 45];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            g = 5.3578547 * X(:, 4).^2 + 0.8356891 * X(:, 2).* X(:, 6) + 37.293239 * X(:, 2) - 40792.141+30665.5386717834 +1;
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1-sqrt(PopObj(:,1))./g);
            % Constraints in objective space
            PopCon(:,1) = max( -(PopObj(:,1).^2 + PopObj(:,2).^2-1), 0);
            % Constraints in decision space
            PopCon(:,2) = + 85.334407 + 0.0056858 * X(:, 3).* X(:, 6) + 0.0006262 * X(:, 2).* X(:, 5) - 0.0022053 * X(:, 4).* X(:, 6) - 92;
            PopCon(:,3) = -85.334407 - 0.0056858 * X(:, 3).* X(:, 6) - 0.0006262 * X(:, 2).* X(:, 5) + 0.0022053 * X(:, 4).* X(:, 6);
            PopCon(:,4) = + 80.51249 + 0.0071317 * X(:, 3).* X(:, 6) + 0.0029955 * X(:, 2).* X(:, 3) + 0.0021813 * X(:, 4).^2 - 110;
            PopCon(:,5) = -80.51249 - 0.0071317 * X(:, 3).* X(:, 6) - 0.0029955 * X(:, 2).* X(:, 3) - 0.0021813 * X(:, 4).^2 + 90;
            PopCon(:,6) = + 9.300961 + 0.0047026 * X(:, 4).* X(:, 6) + 0.0012547 * X(:, 2).* X(:, 4) + 0.0019085 * X(:, 4) .* X(:, 5) - 25;
            PopCon(:,7) = -9.300961 - 0.0047026 * X(:, 4).* X(:, 6) - 0.0012547 * X(:, 2).* X(:, 4) - 0.0019085 * X(:, 4) .* X(:, 5) + 20;
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,2);
            R = R./repmat(sqrt(sum(R.^2,2)),1,2);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end
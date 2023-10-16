classdef DOC6 < PROBLEM
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
            obj.D = 11;
            obj.lower    = [0 -10 * ones(1, 10)];
            obj.upper    = [ 1 10 * ones(1, 10)];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            g_temp = X(:, 2).^2 + X(:, 3).^2 + X(:, 2).* X(:, 3) - 14 * X(:, 2) - 16 * X(:, 3) + (X(:, 4) - 10).^2 + 4 * (X(:, 5) - 5).^2 + ...
            (X(:, 6) - 3).^2 + 2 * (X(:, 7) - 1).^2 + 5 * X(:, 8).^2 + 7 * (X(:, 9) - 11).^2 + 2 * (X(:, 10) - 10).^2 + (X(:, 11) - 7).^2 + 45;
            g = g_temp - 24.3062090681 +1;
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1-sqrt(PopObj(:,1))./g);
            % Constraints in objective space
            c(:,1) = max( -(PopObj(:,1) + PopObj(:,2)-1), 0);
            c(:,2) = max(-(PopObj(:,1) - 0.5).*( PopObj(:,1)+ PopObj(:,2) - 1 - abs(sin(10*pi*(PopObj(:,1) - PopObj(:,2) + 1) ))  ), 0);
            % Constraints in decision space
            c(:,3)  = -105 + 4 * X(:, 2) + 5 * X(:, 3) - 3 * X(:, 8) + 9 * X(:, 9);
            c(:,4)  = 10 * X(:, 2) - 8 * X(:, 3) - 17 * X(:, 8) + 2 * X(:, 9);
            c(:,5)  = -8 * X(:, 2) + 2 * X(:, 3) + 5 * X(:, 10) - 2 * X(:, 11) - 12;
            c(:,6)  = 3 * (X(:, 2) - 2).^2 + 4 * (X(:, 3) - 3).^2 + 2 * X(:, 4).^2 - 7 * X(:, 5) - 120;
            c(:,7)  = 5 * X(:, 2).^2 + 8 * X(:, 3) + (X(:, 4) - 6).^2 - 2 * X(:, 5) - 40;
            c(:,8)  = X(:, 2).^2 + 2 * (X(:, 3) - 2).^2 - 2 * X(:, 2).* X(:, 3) + 14 * X(:, 6) - 6 * X(:, 7);
            c(:,9)  = 0.5 * (X(:, 2) - 8).^2 + 2 * (X(:, 3) - 4).^2 + 3 * X(:, 6).^2 - X(:, 7) - 30;
            c(:,10) = -3 * X(:, 2) + 6 * X(:, 3) + 12 * (X(:, 10) - 8).^2 - 7 * X(:, 11);
            Population  = SOLUTION(X,PopObj,c,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = [(0:0.5/(N-1):0.5),(11:20)/20]';
            R(:,2) = 1 - R(:,1);
        end
    end
end
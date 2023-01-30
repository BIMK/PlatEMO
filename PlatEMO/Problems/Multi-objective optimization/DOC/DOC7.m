classdef DOC7 < PROBLEM
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
            obj.lower    = [0 zeros(1, 10)];
            obj.upper    = [ 1 10 * ones(1, 10)];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X  = varargin{1};
            X  = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            c1 = [-6.089 -17.164 -34.054 -5.914 -24.721 -14.986 -24.1 -10.708 -26.662 -22.179];
            g_temp = sum(X(:,2:11).* (repmat(c1, size(X,1), 1) + log(1E-30 + X(:,2:11)./repmat(1E-30 + sum(X(:,2:11), 2), 1, 10))), 2);
            g = g_temp +47.7648884595 +1;
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1-sqrt(PopObj(:,1))./g);
            % Constraints in objective space
            c(:,1) = max( -(PopObj(:,1) + PopObj(:,2)-1), 0);
            c(:,2) = max(-(PopObj(:,1) - 0.5).*( PopObj(:,1)+ PopObj(:,2) - 1 - abs(sin(10*pi*(PopObj(:,1) - PopObj(:,2) + 1) ))  ), 0);
            c(:,3) = max(- ( abs(- PopObj(:,1) + PopObj(:,2))./sqrt(2) - 0.1./sqrt(2)), 0);
            % Constraints in decision space
            c(:,4) = abs(X(:, 2) + 2 * X(:, 3) + 2 * X(:, 4) + X(:, 7) + X(:, 11) - 2) - 0.0001;
            c(:,5) = abs(X(:, 5) + 2 * X(:, 6) + X(:, 7) + X(:, 8) - 1) - 0.0001;
            c(:,6) = abs(X(:, 4) + X(:, 8) + X(:, 9) + 2 * X(:, 10) + X(:, 11) - 1) - 0.0001;
            Population  = SOLUTION(X,PopObj,c,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = [(0:0.45/(N-1):0.45),(11:20)/20]';
            R(:,2) = 1 - R(:,1);
        end
    end
end
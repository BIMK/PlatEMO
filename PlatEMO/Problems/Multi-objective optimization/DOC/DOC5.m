classdef DOC5 < PROBLEM
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
            obj.D = 8;
            obj.lower    = [0 0 0 0 100 6.3 5.9 4.5];
            obj.upper    = [ 1 1000 40 40 300 6.7 6.4 6.25];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            g_temp = X(:, 2);
            g = g_temp-193.724510070035 +1;
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1-sqrt(PopObj(:,1))./g);
            % Constraints in objective space
            c(:,1) = max( -(PopObj(:,1) + PopObj(:,2)-1), 0);
            c(:,2) = max(- ( PopObj(:,1)+ PopObj(:,2) - 1 - abs(sin(10*pi*(PopObj(:,1) - PopObj(:,2) + 1) ))), 0);
            c(:,3) = max( (PopObj(:,1) - 0.8).*(PopObj(:,2) - 0.6), 0);
            % Constraints in decision space
            c(:,4) = -X(:, 2) + 35 * X(:, 3).^0.6 + 35 * X(:, 4).^0.6;
            c(:,5) = abs(-300 * X(:, 4) + 7500 * X(:, 6) - 7500 * X(:, 7) - 25 * X(:, 5).* X(:, 6) + 25 * X(:, 5).* X(:, 7) + X(:, 4).* X(:, 5)) - 0.0001;
            c(:,6) = abs(100 * X(:, 3) + 155.365 * X(:, 5) + 2500 * X(:, 8) - X(:, 3).* X(:, 5) - 25 * X(:, 5).* X(:, 8) - 15536.5) - 0.0001;
            c(:,7) = abs(-X(:, 6) + log( - X(:, 5) + 900)) - 0.0001;
            c(:,8) = abs(-X(:, 7) + log(X(:, 5) + 300)) - 0.0001;
            c(:,9) = abs(-X(:, 8) + log(-2 * X(:, 5) + 700)) - 0.0001;
            Population = SOLUTION(X,PopObj,c,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = [0:8,16:20]/20;
            R(:,2) = 1 - R(:,1);
        end
    end
end
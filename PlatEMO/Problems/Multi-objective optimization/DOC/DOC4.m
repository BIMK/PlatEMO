classdef DOC4 < PROBLEM
% <multi> <real> <constrained>
% Benchmark MOP with constraints in both decision and objective spaces

%------------------------------- Reference --------------------------------
% Z. Liu and Y. Wang, Handling constrained multiobjective optimization
% problems with constraints in both the decision and objective spaces. IEEE
% Transactions on Evolutionary Computation, 2019, 23(5): 870-884.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
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
            obj.lower    = [0 -10 -10 -10 -10 -10 -10 -10];
            obj.upper    = [ 1 10 10 10 10 10 10 10];
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
            g_temp = (X(:, 2) - 10).^2 + 5 * (X(:, 3) - 12).^2 + X(:, 4).^4 + 3 * (X(:, 5) - 11).^2 + 10 * X(:, 6).^6 + ...
            7 * X(:, 7).^2 + X(:, 8).^4 - 4 * X(:, 7).* X(:, 8) - 10 * X(:, 7) - 8 * X(:, 8);

            g = g_temp-680.6300573745 +1;
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1-sqrt(PopObj(:,1))./g);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            % Constraints in objective space
            c(:,1) = max( -(PopObj(:,1) + PopObj(:,2)-1), 0);
            c(:,2) = max(- ( PopObj(:,1)+ PopObj(:,2) - 1 - abs(sin(10*pi*(PopObj(:,1) - PopObj(:,2) + 1) ))), 0);
            % Constraints in decision space
            c(:, 3) = -127 + 2 * X(:, 2).^2 + 3 * X(:, 3).^4 + X(:, 4) + 4 * X(:, 5).^2 + 5 * X(:, 6);
            c(:, 4) = -282 + 7 * X(:, 2) + 3 * X(:, 3) + 10 * X(:, 4).^2 + X(:, 5) - X(:, 6);
            c(:, 5) = -196 + 23 * X(:, 2) + X(:, 3).^2 + 6 * X(:, 7).^2 - 8 * X(:, 8);
            c(:, 6) = 4 * X(:, 2).^2 + X(:, 3).^2 - 3 * X(:, 2).* X(:, 3) + 2 * X(:, 4).^2 + 5 * X(:, 7) - 11 * X(:, 8);
            PopCon  = c;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = (0:20)/20;
            R(:,2) = 1 - R(:,1);
        end
    end
end
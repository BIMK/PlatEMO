classdef DOC3 < PROBLEM
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
            obj.D = 10;
            obj.lower    = [0 0 0 0 0 0 0 0 0 0.01];
            obj.upper    = [1 1 300 100 200 100 1 100 200 0.03];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            g_temp = -9.* X(:, 6) - 15.* X(:, 9) + 6.* X(:, 2) + 16.* X(:, 3) + 10.* (X(:, 7) + X(:, 8));
            g = (g_temp+400.0551) +1;
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1 - (PopObj(:,1))./g);
            % Constraints in objective space
            c(:,1) = max( -(PopObj(:,1).^2 + PopObj(:,2).^2-1), 0);
            c(:,2) = max(-( abs( (-PopObj(:,1) + PopObj(:,2) -0.5)/sqrt(2)) - 0.1/sqrt(2)), 0);
            c(:,3) = max(-( abs( (-PopObj(:,1) + PopObj(:,2) -0)/sqrt(2)) - 0.1/sqrt(2)), 0);
            c(:,4) = max(-( abs( (-PopObj(:,1) + PopObj(:,2) +0.5)/sqrt(2)) - 0.1/sqrt(2)), 0);

            % Constraints in decision space
            c(:,5)  = X(:, 10).* X(:, 4) + 0.02.* X(:, 7) - 0.025.* X(:, 6);
            c(:,6)  = X(:, 10).* X(:, 5) + 0.02.* X(:, 8) - 0.015.* X(:, 9);
            c(:,7)  = abs(X(:, 2) + X(:, 3) - X(:, 4) - X(:, 5)) - 0.0001;
            c(:,8)  = abs(0.03.* X(:, 2) + 0.01.* X(:, 3) - X(:, 10).* (X(:, 4) + X(:, 5))) - 0.0001;
            c(:,9)  = abs(X(:, 4) + X(:, 7) - X(:, 6)) - 0.0001;
            c(:,10) = abs(X(:, 5) + X(:, 8) - X(:, 9)) - 0.0001;
            Population = SOLUTION(X,PopObj,c,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,2);
            R = R./repmat(sqrt(sum(R.^2,2)),1,2);
            R(0.3403<R(:,1) & R(:,1)<0.4782 | 0.6553<R(:,1) & R(:,1)<0.7553 | 0.8782<R(:,1) & R(:,1)<0.9403,:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = UniformPoint(100,2);
            R = R./repmat(sqrt(sum(R.^2,2)),1,2);
            R(0.3403<R(:,1) & R(:,1)<0.4782 | 0.6553<R(:,1) & R(:,1)<0.7553 | 0.8782<R(:,1) & R(:,1)<0.9403,:) = nan;
        end
    end
end
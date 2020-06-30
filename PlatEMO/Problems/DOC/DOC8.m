classdef DOC8 < PROBLEM
% <problem> <DOC>
% Benchmark MOP with constraints in both decision and objective spaces

%------------------------------- Reference --------------------------------
% Z. Liu and Y. Wang, Handling constrained multiobjective optimization
% problems with constraints in both the decision and objective spaces. IEEE
% Transactions on Evolutionary Computation, 2019.
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
        function obj = DOC8()
            obj.Global.M = 3;
            obj.Global.D = 10;
            obj.Global.lower    = [0 0 500 1000 5000 100 100 100 100 100];
            obj.Global.upper    = [1 1 1000 2000 6000 500 500 500 500 500];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
              % basic multi-objective problem
              g_temp = X(:, 3) + X(:, 4) + X(:, 5);
              g = g_temp-7049.2480205286 +1;
                  
              PopObj(:,1) = (X(:,1).*X(:,2)).*g;
              PopObj(:,2) = (X(:,1).*(1 - X(:,2))).*g;
              PopObj(:,3) = (1-X(:,1)).*g;
            
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            % constraints in objective space
            c(:,1) = max( - (PopObj(:,3) - 0.4).*(PopObj(:,3) - 0.6), 0);
    
            % constraints in decision space
            c(:, 2) = -1 + 0.0025 * (X(:, 6) + X(:, 8));
            c(:, 3) = -1 + 0.0025 * (X(:, 7) + X(:, 9) - X(:, 6));
            c(:, 4) = -1 + 0.01 * (X(:, 10) - X(:, 7));
            c(:, 5) = -X(:, 3).* X(:, 8) + 833.33252 * X(:, 6) + 100 * X(:, 3) - 83333.333;
            c(:, 6) = -X(:, 4).* X(:, 9) + 1250 * X(:, 7) + X(:, 4).* X(:, 6) - 1250 * X(:, 6);
            c(:, 7) = -X(:, 5).* X(:, 10) + 1250000 + X(:, 5).* X(:, 7) - 2500 * X(:, 7);
            PopCon=c;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
             P = UniformPoint(N,3);
             P(0.4<P(:,3)&P(:,3)<0.6,:) = [];
        end
    end
end
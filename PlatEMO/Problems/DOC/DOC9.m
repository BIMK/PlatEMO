classdef DOC9 < PROBLEM
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
        function obj = DOC9()
            obj.Global.M = 3;
            obj.Global.D = 11;
            obj.Global.lower    = [0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1];
            obj.Global.upper    = [1 1 10 10 10 10 10 10 10 10 10];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
              % basic multi-objective problem
              g_temp =  -0.5 * (X(:, 3).* X(:, 6) - X(:, 4).* X(:, 5) + X(:, 5).* X(:, 11) - X(:, 7).* X(:, 11) + X(:, 7).* X(:, 10) - X(:, 8).* X(:, 9));
              g = g_temp +0.8660254038 +1;
    
              PopObj(:,1) = cos(0.5*pi*X(:,1)).*cos(0.5*pi*X(:,2)).*g;
              PopObj(:,2) = cos(0.5*pi*X(:,1)).*sin(0.5*pi*X(:,2)).*g;
              PopObj(:,3) = sin(0.5*pi*X(:,1)).*g;
            
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            % constraints in objective space
            c(:, 1) = max( -( PopObj(:,1).^2 + PopObj(:,2).^2 - 1), 0);
    
            % constraints in decision space
            c(:, 2) = X(:, 5).^2 + X(:, 6).^2 - 1;
            c(:, 3) = X(:, 11).^2 - 1;
            c(:, 4) = X(:, 7).^2 + X(:, 8).^2 - 1;
            c(:, 5) = X(:, 3).^2 + (X(:, 4) - X(:, 11)).^2 - 1;
            c(:, 6) = (X(:, 3) - X(:, 7)).^2 + (X(:, 4) - X(:, 8)).^2 - 1;
            c(:, 7) = (X(:, 3) - X(:, 9)).^2 + (X(:, 4) - X(:, 10)).^2 - 1;
            c(:, 8) = (X(:, 5) - X(:, 7)).^2 + (X(:, 6) - X(:, 8)).^2 - 1;
            c(:, 9) = (X(:, 5) - X(:, 9)).^2 + (X(:, 6) - X(:, 10)).^2 - 1;
            c(:, 10) = X(:, 9).^2 + (X(:, 10) - X(:, 11)).^2 - 1;
            c(:, 11) = X(:, 4).* X(:, 5) - X(:, 3).* X(:, 6);
            c(:, 12) = -X(:, 5).* X(:, 11);
            c(:, 13) = X(:, 7).* X(:, 11);
            c(:, 14) = X(:, 8).* X(:, 9) - X(:, 7).* X(:, 10);
            PopCon=c;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,2);
            P = P./repmat(sqrt(sum(P.^2,2)),1,2);
            P = [P,zeros(size(P,1),1)];
        end
    end
end
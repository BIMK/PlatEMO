classdef DOC1 < PROBLEM
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
        function obj = DOC1()
            obj.Global.M = 2;
            obj.Global.D = 6;
            obj.Global.lower    = [0 78 33 27 27 27];
            obj.Global.upper    = [1 102 45 45 45 45];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
             g      = 5.3578547 * X(:, 4).^2 + 0.8356891 * X(:, 2).* X(:, 6) + 37.293239 * X(:, 2) - 40792.141+30665.5386717834 +1;
             PopObj(:,1) = X(:,1);
             PopObj(:,2) = g.*(1-sqrt(PopObj(:,1))./g);
            
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);

            % constraints in objective space
           PopCon(:,1) = max( -(PopObj(:,1).^2 + PopObj(:,2).^2-1), 0);
    
            % constraints in decision space
            PopCon(:, 2) = + 85.334407 + 0.0056858 * X(:, 3).* X(:, 6) + 0.0006262 * X(:, 2).* X(:, 5) - 0.0022053 * X(:, 4).* X(:, 6) - 92;
            PopCon(:, 3) = -85.334407 - 0.0056858 * X(:, 3).* X(:, 6) - 0.0006262 * X(:, 2).* X(:, 5) + 0.0022053 * X(:, 4).* X(:, 6);
            PopCon(:, 4) = + 80.51249 + 0.0071317 * X(:, 3).* X(:, 6) + 0.0029955 * X(:, 2).* X(:, 3) + 0.0021813 * X(:, 4).^2 - 110;
            PopCon(:, 5) = -80.51249 - 0.0071317 * X(:, 3).* X(:, 6) - 0.0029955 * X(:, 2).* X(:, 3) - 0.0021813 * X(:, 4).^2 + 90;
            PopCon(:, 6) = + 9.300961 + 0.0047026 * X(:, 4).* X(:, 6) + 0.0012547 * X(:, 2).* X(:, 4) + 0.0019085 * X(:, 4) .* X(:, 5) - 25;
            PopCon(:, 7) = -9.300961 - 0.0047026 * X(:, 4).* X(:, 6) - 0.0012547 * X(:, 2).* X(:, 4) - 0.0019085 * X(:, 4) .* X(:, 5) + 20;
            
         end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,2);
            P = P./repmat(sqrt(sum(P.^2,2)),1,2);
        end
    end
end
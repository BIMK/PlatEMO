classdef DOC2 < PROBLEM
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
            obj.D = 16;
            obj.lower    = [0 zeros(1, 15)];
            obj.upper    = [1 10 * ones(1, 15)];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            
            [popsize,~] = size(X);

            b  = [-40 -2 -0.25 -4 -4 -1 -40 -60 5 1];
            c1 = [30 -20 -10 32 -10;
                 -20 39 -6 -31 32;
                 -10 -6 10 -6 -10;
                  32 -31 -6 39 -20;
                 -10 32 -10 -20 30];
            d  = [4 8 10 6 2];

            g_temp = sum(repmat(c1(1:5, 1)', popsize, 1).* X(:, 12:16), 2).* X(:, 12) + sum(repmat(c1(1:5, 2)', popsize, 1).* X(:, 12:16), 2).* X(:, 13)...
                     + sum(repmat(c1(1:5, 3)', popsize, 1).* X(:, 12:16), 2).* X(:, 14) + sum(repmat(c1(1:5, 4)', popsize, 1).* X(:, 12:16), 2).* X(:, 15)...
                     + sum(repmat(c1(1:5, 5)', popsize, 1).* X(:, 12:16), 2).* X(:, 16) + 2 * sum(repmat(d, popsize, 1).* X(:, 12:16).^3, 2)...
                     - sum(repmat(b, popsize, 1).* X(:, 2:11), 2);
            g = (g_temp-32.6555929502) +1;

            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1- (PopObj(:,1)).^(1/3)./g);
            % Constraints in objective space
            PopCon(:,1) = max( -(sqrt(PopObj(:,1)) + PopObj(:,2)-1), 0);
            d1(:,1) = max(((PopObj(:,1)-1/8).^2 + (PopObj(:,2) -1+sqrt(1/8) ).^2 - 0.15*0.15),0);
            d1(:,2) = max( ((PopObj(:,1)-1/2).^2 + (PopObj(:,2) -1+sqrt(1/2)).^2 - 0.15*0.15),0);  
            d1(:,3) = max(((PopObj(:,1)-7/8).^2 + (PopObj(:,2) - 1+sqrt(7/8)).^2 - 0.15*0.15),0);    
            PopCon(:,2) = min(d1, [], 2);
            
            a = [-16 2 0 1 0;
            0 -2 0 0.4 2;
            -3.5 0 2 0 0;
            0 -2 0 -4 -1;
            0 -9 -2 1 -2.8;
            2 0 -4 0 0;
            -1 -1 -1 -1 -1;
            -1 -2 -3 -2 -1;
            1 2 3 4 5;
            1 1 1 1 1];
            
            c1 = [30 -20 -10 32 -10;
            -20 39 -6 -31 32;
            -10 -6 10 -6 -10;
            32 -31 -6 39 -20;
            -10 32 -10 -20 30];
            d = [4 8 10 6 2];
            e = [-15 -27 -36 -18 -12];
            % Constraints in decision space
            PopCon(:,3) = -2 * sum(repmat(c1(1:5, 1)', popsize, 1).* X(:, 12:16), 2) - 3 * d(1).* X(:, 12).^2 - e(1) + sum(repmat(a(1:10, 1)', popsize, 1).* X(:, 2:11), 2);
            PopCon(:,4) = -2 * sum(repmat(c1(1:5, 2)', popsize, 1).* X(:, 12:16), 2) - 3 * d(2).* X(:, 13).^2 - e(2) + sum(repmat(a(1:10, 2)', popsize, 1).* X(:, 2:11), 2);
            PopCon(:,5) = -2 * sum(repmat(c1(1:5, 3)', popsize, 1).* X(:, 12:16), 2) - 3 * d(3).* X(:, 14).^2 - e(3) + sum(repmat(a(1:10, 3)', popsize, 1).* X(:, 2:11), 2);
            PopCon(:,6) = -2 * sum(repmat(c1(1:5, 4)', popsize, 1).* X(:, 12:16), 2) - 3 * d(4).* X(:, 15).^2 - e(4) + sum(repmat(a(1:10, 4)', popsize, 1).* X(:, 2:11), 2);
            PopCon(:,7) = -2 * sum(repmat(c1(1:5, 5)', popsize, 1).* X(:, 12:16), 2) - 3 * d(5).* X(:, 16).^2 - e(5) + sum(repmat(a(1:10, 5)', popsize, 1).* X(:, 2:11), 2);
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - sqrt(R(:,1));
            R(R(:,1)<0.05 | 0.2202<R(:,1) & R(:,1)<0.3830 | 0.6247<R(:,1) & R(:,1)<0.7440,:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R(:,1) = linspace(0,1,100)';
            R(:,2) = 1 - sqrt(R(:,1));
            R(R(:,1)<0.05 | 0.2202<R(:,1) & R(:,1)<0.3830 | 0.6247<R(:,1) & R(:,1)<0.7440,:) = nan;
        end
    end
end
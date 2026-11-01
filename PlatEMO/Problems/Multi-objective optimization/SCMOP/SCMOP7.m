classdef SCMOP7 < PROBLEM
% <2026> <multi/many> <real> <large/none> <expensive/none> <sparse> <constrained>
% Benchmark CMOP with sparse Pareto optimal solutions
% theta_prime --- 0.1 --- Sparsity level
% Tns         --- 0.2 --- Difficulty parameter of nonsparse constraints
% Ts          --- 0.2 --- Difficulty parameter of sparse constraints

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, Y. Zhang, and X. Zhang. An evolutionary algorithm with
% adaptive search space switching and constraint-guided operators for
% sparse constrained multiobjective optimization. IEEE Transactions on
% Systems, Man, and Cybernetics: Systems, 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        theta = 0.1;    % theta_prime Sparsity level  
        C1    = 0.2;    % Tns Difficulty parameter of nonsparse constraints
        C2    = 0.3;    % Ts Difficulty parameter of sparse constraints
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.theta, obj.C1, obj.C2] = obj.ParameterSet(0.1,0.2,0.3);
            if isempty(obj.M)
                obj.M = 2;
            end
            if isempty(obj.D)
                obj.D = 100;
            end
            obj.lower    = [zeros(1,obj.M-1)+0,zeros(1,obj.D-obj.M+1)-1];
            obj.upper    = [zeros(1,obj.M-1)+1,zeros(1,obj.D-obj.M+1)+2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            g      = sum(g2(X(:,obj.M:end),optx1(X(:,obj.M:end))),2);
            PopObj = repmat(1+g/(obj.D-obj.M+1),1,obj.M).*fliplr(cumprod([ones(size(X,1),1),cos(X(:,1:obj.M-1)*pi/2)],2)).*[ones(size(X,1),1),sin(X(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            K           = ceil(obj.theta*(obj.D-obj.M+1));
            PopCon(:,1) = mean(g3(X(:,obj.M:obj.M+K-1),optx2(X(:,obj.M:obj.M+K-1))),2) - obj.C1;
            PopCon(:,2) = mean(g1(X(:,obj.M+K:end),0),2) - obj.C2;    
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {cos(a)*cos(a'),cos(a)*sin(a'),sin(a)*ones(size(a'))};
            else
                R = [];
            end
        end
    end
end

function x = optx1(x)
    x = min(0.5,floor(abs(x) * 10) / 10);
end

function x = optx2(x)
    x = max(0.5, floor(abs(x) * 10) / 10);
end

function g = g1(x,t)
    g = (x-t).^2;
end

function g = g2(x,t)
    g = 2*(x-t).^2 + sin(2*pi*(x-t)).^2;
end

function g = g3(x,t)
    g = 4-(x-t)-4./exp(100*(x-t).^2);
end
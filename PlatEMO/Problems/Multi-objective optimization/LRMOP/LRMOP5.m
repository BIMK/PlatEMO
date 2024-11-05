classdef LRMOP5 < PROBLEM
% <2024> <multi/many> <real> <large/none> <robust> <sparse/none>
% Sparse robust multi-objective optimization problem
% theta --- 0.1 --- Sparsity of the Pareto set
% H     ---  50 --- Number of disturbances

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, L. Zhang, K. C. Tan, and X. Zhang, An evolutionary
% algorithm for solving large-scale robust multi-objective optimization
% problems, IEEE Transactions on Evolutionary Computation, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        theta = 0.1;    % Sparsity of the Pareto set
        H     = 50;     % Number of disturbances
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.theta,obj.H] = obj.ParameterSet(0.1,50);
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 100; end
            obj.lower    = [zeros(1,obj.M-1)+0,zeros(1,obj.D-obj.M+1)-1];
            obj.upper    = [zeros(1,obj.M-1)+1,zeros(1,obj.D-obj.M+1)+2];
            obj.encoding = ones(1,obj.D);
        end     
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [~,D] = size(X);
            M = obj.M;
            k = ceil(obj.theta*(obj.D-obj.M+1));
            K = k;
            a= random('Normal',1,0.2);
            g = sum(g1(X(:,2:end),a).*g2(X(:,2:end),0),2) + abs(K-sum(X(:,2:end)~=0,2));
            PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(size(X,1),1),cos(X(:,1:obj.M-1)*pi/2)],2)).*[ones(size(X,1),1),sin(X(:,obj.M-1:-1:1)*pi/2)];
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
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            switch metName
                case {'Mean_IGD','Mean_HV','Worst_IGD','Worst_HV'}
                    score = feval(metName,Population,obj);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
        %% Perturb solutions multiple times
        function PopX = Perturb(obj,PopDec,N)
            if nargin < 3; N = obj.H; end
            PopX = [];
            for i = 1 : N
                PopX = [PopX,SOLUTION(PopDec,obj.CalObj(PopDec),obj.CalCon(PopDec))];
            end
        end
    end
end
function g = g2(x,t)   
     g = 4-(x-t)-4./exp(100*(x-t).^2);
end

function g = g1(x,t)  
   g = (3+x-4*t*t*t).^2;
end
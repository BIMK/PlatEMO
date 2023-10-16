classdef MOEADM2M_F1 < PROBLEM
% <multi> <real> <large/none>
% Benchmark MOP for testing MOEA/D-M2M

%------------------------------- Reference --------------------------------
% H. Liu, F. Gu, and Q. Zhang, Decomposition of a multiobjective
% optimization problem into a number of simple multiobjective subproblems,
% IEEE Transactions on Evolutionary Computation, 2014, 18(3): 450-455.
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
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            t = X(:,2:end) - repmat(sin(pi/2*X(:,1)),1,size(X,2)-1);
            g = sin(pi*X(:,1)).*sum(-0.9*t.^2+abs(t).^0.6,2);
            PopObj(:,1) = (1+g).*X(:,1);
            PopObj(:,2) = (1+g).*(1-sqrt(X(:,1)));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - sqrt(R(:,1));
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end
classdef IMMOEA_F9 < PROBLEM
% <multi> <real> <large/none>
% Benchmark MOP for testing IM-MOEA

%------------------------------- Reference --------------------------------
% R. Cheng, Y. Jin, K. Narukawa, and B. Sendhoff, A multiobjective
% evolutionary algorithm using Gaussian process-based inverse modeling,
% IEEE Transactions on Evolutionary Computation, 2015, 19(6): 838-856.
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
            if isempty(obj.D); obj.D = 30; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = [1,zeros(1,obj.D-1)+10];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            t = X(:,2:obj.D).^(1./(1+3*repmat(2:obj.D,size(X,1),1)/obj.D)) - repmat(X(:,1),1,obj.D-1);
            g = sum(t.^2/4000,2) - prod(cos(t./repmat(sqrt(1:obj.D-1),size(X,1),1)),2) + 2;
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1-sqrt(PopObj(:,1)./g));
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
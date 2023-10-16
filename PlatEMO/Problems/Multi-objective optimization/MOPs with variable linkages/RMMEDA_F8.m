classdef RMMEDA_F8 < PROBLEM
% <multi> <real> <large/none>
% Benchmark MOP for testing RM-MEDA

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, and Y. Jin, RM-MEDA: A regularity model-based
% multiobjective estimation of distribution algorithm, IEEE Transactions on
% Evolutionary Computation, 2008, 12(1): 41-63.
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
            obj.M = 3;
            if isempty(obj.D); obj.D = 30; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            g = sum((X(:,3:end).^2-repmat(X(:,1),1,size(X,2)-2)).^2,2);
            PopObj(:,1) = cos(pi/2*X(:,1)).*cos(pi/2*X(:,2)).*(1+g);
            PopObj(:,2) = cos(pi/2*X(:,1)).*sin(pi/2*X(:,2)).*(1+g);
            PopObj(:,3) = sin(pi/2*X(:,1)).*(1+g);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,3);
            R = R./repmat(sqrt(sum(R.^2,2)),1,3);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            a = linspace(0,pi/2,10)';
            R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
        end
    end
end
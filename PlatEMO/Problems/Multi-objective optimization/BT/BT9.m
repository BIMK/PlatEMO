classdef BT9 < PROBLEM
% <multi> <real> <large/none>
% Benchmark MOP with bias feature

%------------------------------- Reference --------------------------------
% H. Li, Q. Zhang, and J. Deng, Biased multiobjective optimization and
% decomposition algorithm, IEEE Transactions on Cybernetics, 2017, 47(1):
% 52-66.
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
            [N,D] = size(X);
            I1    = 3 : 3 : D;
            I2    = 4 : 3 : D;
            I3    = 5 : 3 : D;
            Y     = X - sin(repmat(1:D,N,1)*pi/2/D);
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + sum(Y(:,I1).^2+(1-exp(-Y(:,I1).^2/1e-9))/5,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + sum(Y(:,I2).^2+(1-exp(-Y(:,I2).^2/1e-9))/5,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + sum(Y(:,I3).^2+(1-exp(-Y(:,I3).^2/1e-9))/5,2);
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
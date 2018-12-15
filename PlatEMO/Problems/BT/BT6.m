classdef BT6 < PROBLEM
% <problem> <BT>
% Benchmark MOP with bias feature

%------------------------------- Reference --------------------------------
% H. Li, Q. Zhang, and J. Deng, Biased multiobjective optimization and
% decomposition algorithm, IEEE Transactions on Cybernetics, 2017, 47(1):
% 52-66.
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
        function obj = BT6()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 30;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,D] = size(X);
            I1    = 2 : 2 : D;
            I2    = 3 : 2 : D;
            Y     = X - repmat(X(:,1),1,D).^(0.5+1.5*repmat(0:D-1,N,1)/(D-1));
            PopObj(:,1) = X(:,1)         + sum(Y(:,I1).^2+(1-exp(-Y(:,I1).^2/1e-4))/5,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + sum(Y(:,I2).^2+(1-exp(-Y(:,I2).^2/1e-4))/5,2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1).^0.5;
        end
    end
end
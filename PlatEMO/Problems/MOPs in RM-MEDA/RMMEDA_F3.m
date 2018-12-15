classdef RMMEDA_F3 < PROBLEM
% <problem> <RMMEDA>
% Benchmark MOP for RM-MEDA

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, and Y. Jin, RM-MEDA: A regularity model-based
% multiobjective estimation of distribution algorithm, IEEE Transactions on
% Evolutionary Computation, 2008, 12(1): 41-63.
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
        function obj = RMMEDA_F3()
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
            g = 1 + 9*(sum((X(:,2:end)-repmat(X(:,1),1,size(X,2)-1)).^2,2)/9).^0.25;
            PopObj(:,1) = 1 - exp(-4*X(:,1)).*sin(6*pi*X(:,1)).^6;
            PopObj(:,2) = g.*(1-(PopObj(:,1)./g).^2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            minf1  = min(1-exp(-4*(0:1e-6:1)).*(sin(6*pi*(0:1e-6:1))).^6);
            P(:,1) = (minf1:(1-minf1)/(N-1):1)';
            P(:,2) = 1 - P(:,1).^2;
        end
    end
end
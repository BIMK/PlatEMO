classdef MOEADDE_F3 < PROBLEM
% <problem> <MOEADDE>
% Benchmark MOP for MOEA/D-DE

%------------------------------- Reference --------------------------------
% H. Li and Q. Zhang, Multiobjective optimization problems with complicated
% Pareto sets, MOEA/D and NSGA-II, IEEE Transactions on Evolutionary
% Computation, 2009, 13(2): 284-302.
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
        function obj = MOEADDE_F3()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 30;
            end
            obj.Global.lower    = [0,-ones(1,obj.Global.D-1)];
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,D] = size(X);
            J1    = 3 : 2 : D;
            J2    = 2 : 2 : D;
            PopObj(:,1) = X(:,1)         + 2*mean((X(:,J1)-0.8*repmat(X(:,1),1,length(J1)).*cos(repmat(6*pi*X(:,1),1,length(J1))+repmat(J1*pi/D,N,1))).^2,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean((X(:,J2)-0.8*repmat(X(:,1),1,length(J2)).*sin(repmat(6*pi*X(:,1),1,length(J2))+repmat(J2*pi/D,N,1))).^2,2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1).^0.5;
        end
    end
end
classdef MOEADDE_F6 < PROBLEM
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
        function obj = MOEADDE_F6()
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = [0,0,zeros(1,obj.Global.D-2)-2];
            obj.Global.upper    = [1,1,zeros(1,obj.Global.D-2)+2];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,D] = size(X);
            J1    = 4 : 3 : D;
            J2    = 5 : 3 : D;
            J3    = 3 : 3 : D;
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean((X(:,J1)-2*repmat(X(:,2),1,length(J1)).*sin(repmat(2*pi*X(:,1),1,length(J1))+repmat(J1*pi/D,N,1))).^2,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean((X(:,J2)-2*repmat(X(:,2),1,length(J2)).*sin(repmat(2*pi*X(:,1),1,length(J2))+repmat(J2*pi/D,N,1))).^2,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean((X(:,J3)-2*repmat(X(:,2),1,length(J3)).*sin(repmat(2*pi*X(:,1),1,length(J3))+repmat(J3*pi/D,N,1))).^2,2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,3);
            P = P./repmat(sqrt(sum(P.^2,2)),1,3);
        end
    end
end
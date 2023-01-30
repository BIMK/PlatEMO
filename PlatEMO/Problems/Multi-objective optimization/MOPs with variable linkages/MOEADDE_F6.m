classdef MOEADDE_F6 < PROBLEM
% <multi> <real> <large/none>
% Benchmark MOP for testing MOEA/D-DE

%------------------------------- Reference --------------------------------
% H. Li and Q. Zhang, Multiobjective optimization problems with complicated
% Pareto sets, MOEA/D and NSGA-II, IEEE Transactions on Evolutionary
% Computation, 2009, 13(2): 284-302.
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
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = [0,0,zeros(1,obj.D-2)-2];
            obj.upper    = [1,1,zeros(1,obj.D-2)+2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            J1 = 4 : 3 : obj.D;
            J2 = 5 : 3 : obj.D;
            J3 = 3 : 3 : obj.D;
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean((X(:,J1)-2*repmat(X(:,2),1,length(J1)).*sin(repmat(2*pi*X(:,1),1,length(J1))+repmat(J1*pi/obj.D,size(X,1),1))).^2,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean((X(:,J2)-2*repmat(X(:,2),1,length(J2)).*sin(repmat(2*pi*X(:,1),1,length(J2))+repmat(J2*pi/obj.D,size(X,1),1))).^2,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean((X(:,J3)-2*repmat(X(:,2),1,length(J3)).*sin(repmat(2*pi*X(:,1),1,length(J3))+repmat(J3*pi/obj.D,size(X,1),1))).^2,2);
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
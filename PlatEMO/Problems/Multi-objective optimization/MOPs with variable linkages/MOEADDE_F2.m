classdef MOEADDE_F2 < PROBLEM
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
            obj.M = 2;
            if isempty(obj.D); obj.D = 30; end
            obj.lower    = [0,-ones(1,obj.D-1)];
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            J1 = 3 : 2 : obj.D;
            J2 = 2 : 2 : obj.D;
            PopObj(:,1) = X(:,1)         + 2*mean((X(:,J1)-sin(repmat(6*pi*X(:,1),1,length(J1))+repmat(J1*pi/obj.D,size(X,1),1))).^2,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean((X(:,J2)-sin(repmat(6*pi*X(:,1),1,length(J2))+repmat(J2*pi/obj.D,size(X,1),1))).^2,2);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - R(:,1).^0.5;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end
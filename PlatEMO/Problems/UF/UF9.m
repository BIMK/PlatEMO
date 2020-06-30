classdef UF9 < PROBLEM
% <problem> <UF>
% Unconstrained benchmark MOP

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, S. Zhao, P. N. Suganthan, W. Liu, and S. Tiwari,
% Multiobjective optimization test instances for the CEC 2009 special
% session and competition, School of CS & EE, University of Essex, Working
% Report CES-487, 2009.
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
        function obj = UF9()
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 30;
            end
            obj.Global.lower    = [0,0,zeros(1,obj.Global.D-2)-2];
            obj.Global.upper    = [1,1,zeros(1,obj.Global.D-2)+2];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            D  = size(X,2);
            J1 = 4 : 3 : D;
            J2 = 5 : 3 : D;
            J3 = 3 : 3 : D;
            Y  = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
            PopObj(:,1) = 0.5*(max(0,1.1*(1-4*(2*X(:,1)-1).^2))+2*X(:,1)).*X(:,2)   + 2*mean(Y(:,J1).^2,2);
            PopObj(:,2) = 0.5*(max(0,1.1*(1-4*(2*X(:,1)-1).^2))-2*X(:,1)+2).*X(:,2) + 2*mean(Y(:,J2).^2,2);
            PopObj(:,3) = 1-X(:,2)                                                  + 2*mean(Y(:,J3).^2,2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,3);
            P(P(:,1)>(1-P(:,3))/4 & P(:,1)<(1-P(:,3))*3/4,:) = [];
        end
    end
end
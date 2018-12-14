classdef CF9 < PROBLEM
% <problem> <CF>
% Constrained benchmark MOP

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
        function obj = CF9()
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
            D  = size(X,2);
            J1 = 4 : 3 : D;
            J2 = 5 : 3 : D;
            J3 = 3 : 3 : D;
            Y  = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(Y(:,J1).^2,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(Y(:,J2).^2,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(Y(:,J3).^2,2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            PopCon = 1 - (PopObj(:,1).^2+PopObj(:,2).^2)./(1-PopObj(:,3).^2) + 3*sin(2*pi*((PopObj(:,1).^2-PopObj(:,2).^2)./(1-PopObj(:,3).^2)+1));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,3);
            P = P./repmat(sqrt(sum(P.^2,2)),1,3);
            P(1e-5<P(:,1) & P(:,1)<sqrt((1-P(:,3).^2)/4) | sqrt((1-P(:,3).^2)/2)<P(:,1) & P(:,1)<sqrt(3*(1-P(:,3).^2)/4),:) = [];
        end
    end
end
classdef CF4 < PROBLEM
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
        function obj = CF4()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = [0,zeros(1,obj.Global.D-1)-2];
            obj.Global.upper    = [1,zeros(1,obj.Global.D-1)+2];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            D  = size(X,2);
            J1 = 3 : 2 : D;
            J2 = 2 : 2 : D;
            Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
            h           = Y.^2;
            temp        = Y(:,2) < 3/2*(1-sqrt(1/2));
            h(temp,2)   = abs(Y(temp,2));
            h(~temp,2)  = 0.125 + (Y(~temp,2)-1).^2;
            PopObj(:,1) = X(:,1)   + sum(h(:,J1),2);
            PopObj(:,2) = 1-X(:,1) + sum(h(:,J2),2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            t      = X(:,2) - sin(6*pi*X(:,1)+2*pi/size(X,2)) - 0.5*X(:,1) + 0.25;
            PopCon = -t./(1+exp(4*abs(t)));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1);
            temp1  = 0.5<P(:,1) & P(:,1)<=0.75;
            temp2  = 0.75<P(:,1);
            P(temp1,2) = -0.5*P(temp1,1)+3/4;
            P(temp2,2) = 1 - P(temp2,1)+0.125;
        end
    end
end
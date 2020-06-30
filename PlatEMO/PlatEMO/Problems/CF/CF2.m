classdef CF2 < PROBLEM
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
        function obj = CF2()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = [0,zeros(1,obj.Global.D-1)-1];
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            D  = size(X,2);
            J1 = 3 : 2 : D;
            J2 = 2 : 2 : D;
            PopObj(:,1) = X(:,1)         + 2*mean((X(:,J1)-sin(6*pi*repmat(X(:,1),1,length(J1))+repmat(J1,size(X,1),1)*pi/D)).^2,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean((X(:,J2)-cos(6*pi*repmat(X(:,1),1,length(J2))+repmat(J2,size(X,1),1)*pi/D)).^2,2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            t      = PopObj(:,2) + sqrt(PopObj(:,1)) - sin(2*pi*(sqrt(PopObj(:,1))-PopObj(:,2)+1)) - 1;
            PopCon = -t./(1+exp(4*abs(t)));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - sqrt(P(:,1));
            P(0<P(:,1) & P(:,1)<1/16 | 1/4<P(:,1) & P(:,1)<9/16,:) = [];
        end
    end
end
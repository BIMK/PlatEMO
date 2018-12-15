classdef UF2 < PROBLEM
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
        function obj = UF2()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 30;
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
            Y       = zeros(size(X));
            X1      = repmat(X(:,1),1,length(J1));
            Y(:,J1) = X(:,J1)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J1,size(X,1),1)*pi/D)+0.6*X1).*cos(6*pi*X1+repmat(J1,size(X,1),1)*pi/D);
            X1      = repmat(X(:,1),1,length(J2));
            Y(:,J2) = X(:,J2)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J2,size(X,1),1)*pi/D)+0.6*X1).*sin(6*pi*X1+repmat(J2,size(X,1),1)*pi/D); 
            PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1).^0.5;
        end
    end
end
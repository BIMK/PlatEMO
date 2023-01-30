classdef CF7 < PROBLEM
% <multi> <real> <large/none> <constrained>
% Constrained benchmark MOP

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, S. Zhao, P. N. Suganthan, W. Liu, and S. Tiwari,
% Multiobjective optimization test instances for the CEC 2009 special
% session and competition, School of CS & EE, University of Essex, Working
% Report CES-487, 2009.
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
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = [0,zeros(1,obj.D-1)-2];
            obj.upper    = [1,zeros(1,obj.D-1)+2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            D  = size(X,2);
            J1 = 3 : 2 : D;
            J2 = 2 : 2 : D;
            Y           = zeros(size(X));
            Y(:,J1)     = X(:,J1) - cos(6*pi*repmat(X(:,1),1,length(J1))+repmat(J1,size(X,1),1)*pi/D);
            Y(:,J2)     = X(:,J2) - sin(6*pi*repmat(X(:,1),1,length(J2))+repmat(J2,size(X,1),1)*pi/D);
            h           = 2*Y.^2 - cos(4*pi*Y) + 1;
            h(:,[2,4])  = Y(:,[2,4]).^2;
            PopObj(:,1) = X(:,1)        + sum(h(:,J1),2);
            PopObj(:,2) = (1-X(:,1)).^2 + sum(h(:,J2),2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopCon(:,1) = -X(:,2) + sin(6*pi*X(:,1)+2*pi/size(X,2)) + sign(0.5*(1-X(:,1))-(1-X(:,1)).^2).*sqrt(abs(0.5*(1-X(:,1))-(1-X(:,1)).^2));
            PopCon(:,2) = -X(:,4) + sin(6*pi*X(:,1)+4*pi/size(X,2)) + sign(0.25*sqrt(1-X(:,1))-0.5*(1-X(:,1))).*sqrt(abs(0.25*sqrt(1-X(:,1))-0.5*(1-X(:,1))));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = (1-R(:,1)).^2;
            temp1  = 0.5<R(:,1) & R(:,1)<=0.75;
            temp2  = 0.75<R(:,1);
            R(temp1,2) = 0.5*(1-R(temp1,1));
            R(temp2,2) = 0.25*sqrt(1-R(temp2,1));
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end
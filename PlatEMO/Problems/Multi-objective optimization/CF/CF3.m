classdef CF3 < PROBLEM
% <multi> <real> <large/none> <constrained>
% Constrained benchmark MOP

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, S. Zhao, P. N. Suganthan, W. Liu, and S. Tiwari,
% Multiobjective optimization test instances for the CEC 2009 special
% session and competition, School of CS & EE, University of Essex, Working
% Report CES-487, 2009.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
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
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            D  = size(X,2);
            J1 = 3 : 2 : D;
            J2 = 2 : 2 : D;
            Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
            PopObj(:,1) = X(:,1)      + 2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
            PopObj(:,2) = 1-X(:,1).^2 + 2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            PopCon = 1 - PopObj(:,2) - PopObj(:,1).^2 + sin(2*pi*(PopObj(:,1).^2-PopObj(:,2)+1));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - R(:,1).^2;
            R(0<R(:,1) & R(:,1)<1/2 | sqrt(1/2)<R(:,1) & R(:,1)<sqrt(3/4),:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R(:,1) = linspace(0,1,100)';
            R(:,2) = 1 - R(:,1).^2;
            R(0<R(:,1) & R(:,1)<1/2 | sqrt(1/2)<R(:,1) & R(:,1)<sqrt(3/4),:) = nan;
        end
    end
end
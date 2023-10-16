classdef CF10 < PROBLEM
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
            obj.M = 3;
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = [0,0,zeros(1,obj.D-2)-2];
            obj.upper    = [1,1,zeros(1,obj.D-2)+2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X  = varargin{1};
            X  = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            D  = size(X,2);
            J1 = 4 : 3 : D;
            J2 = 5 : 3 : D;
            J3 = 3 : 3 : D;
            Y  = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(4*Y(:,J1).^2-cos(8*pi*Y(:,J1))+1,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(4*Y(:,J2).^2-cos(8*pi*Y(:,J2))+1,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(4*Y(:,J3).^2-cos(8*pi*Y(:,J3))+1,2);
            PopCon      = 1 - (PopObj(:,1).^2+PopObj(:,2).^2)./(1-PopObj(:,3).^2) + sin(2*pi*((PopObj(:,1).^2-PopObj(:,2).^2)./(1-PopObj(:,3).^2)+1));
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,3);
            R = R./repmat(sqrt(sum(R.^2,2)),1,3);
            R(1e-5<R(:,1) & R(:,1)<sqrt((1-R(:,3).^2)/4) | sqrt((1-R(:,3).^2)/2)<R(:,1) & R(:,1)<sqrt(3*(1-R(:,3).^2)/4),:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            a = linspace(0,pi/2,20)';
            x = sin(a)*cos(a');
            y = sin(a)*sin(a');
            z = cos(a)*ones(size(a'));
            infes    = 1e-5<x & x<sqrt((1-z.^2)/4) | sqrt((1-z.^2)/2)<x & x<sqrt(3*(1-z.^2)/4);
            z(infes) = nan;
            R = {x,y,z};
        end
    end
end
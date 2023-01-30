classdef CF2 < PROBLEM
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
            obj.lower    = [0,zeros(1,obj.D-1)-1];
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X  = varargin{1};
            X  = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            D  = size(X,2);
            J1 = 3 : 2 : D;
            J2 = 2 : 2 : D;
            PopObj(:,1) = X(:,1)         + 2*mean((X(:,J1)-sin(6*pi*repmat(X(:,1),1,length(J1))+repmat(J1,size(X,1),1)*pi/D)).^2,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean((X(:,J2)-cos(6*pi*repmat(X(:,1),1,length(J2))+repmat(J2,size(X,1),1)*pi/D)).^2,2);
            t           = PopObj(:,2) + sqrt(PopObj(:,1)) - sin(2*pi*(sqrt(PopObj(:,1))-PopObj(:,2)+1)) - 1;
            PopCon      = -t./(1+exp(4*abs(t)));
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - sqrt(R(:,1));
            R(0<R(:,1) & R(:,1)<1/16 | 1/4<R(:,1) & R(:,1)<9/16,:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R(:,1) = linspace(0,1,100)';
            R(:,2) = 1 - sqrt(R(:,1));
            R(0<R(:,1) & R(:,1)<1/16 | 1/4<R(:,1) & R(:,1)<9/16,:) = nan;
        end
    end
end
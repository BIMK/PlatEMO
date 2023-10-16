classdef CF1 < PROBLEM
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
            obj.lower    = zeros(1,obj.D);
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
            PopObj(:,1) = X(:,1)   + 2*mean((X(:,J1)-repmat(X(:,1),1,length(J1)).^(0.5*(1+3*(repmat(J1,size(X,1),1)-2)/(D-2)))).^2,2);
            PopObj(:,2) = 1-X(:,1) + 2*mean((X(:,J2)-repmat(X(:,1),1,length(J2)).^(0.5*(1+3*(repmat(J2,size(X,1),1)-2)/(D-2)))).^2,2);
            PopCon      = 1 - PopObj(:,1) - PopObj(:,2) + abs(sin(10*pi*(PopObj(:,1)-PopObj(:,2)+1)));
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = (0:1/20:1)';
            R(:,2) = 1 - R(:,1);
        end
    end
end
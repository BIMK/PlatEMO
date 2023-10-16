classdef LIRCMOP13 < PROBLEM
% <multi> <real> <large/none> <constrained>
% Constrained benchmark MOP with large infeasible regions

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Huang, Y. Fang, Y. You, J. Mo, C. Wei, and E.
% Goodman, An improved epsilon constraint-handling method in MOEA/D for
% CMOPs with large infeasible regions, Soft Computing, 2019, 23:
% 12491-12510.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenji Li
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 3;
            if isempty(obj.D); obj.D = 30; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            [popsize,variable_length] = size(X);
            sum1 = zeros(popsize,1);
            for j = 3 : variable_length
                sum1 = sum1+10*(X(:,j)-0.5).^2;
            end
            PopObj(:,1) = (1.7057+sum1).*cos(0.5*pi*X(:,1)).*cos(0.5*pi*X(:,2));
            PopObj(:,2) = (1.7057+sum1).*cos(0.5*pi*X(:,1)).*sin(0.5*pi*X(:,2));
            PopObj(:,3) = (1.7057+sum1).*sin(0.5*pi*X(:,1));
            gx          =  PopObj(:,1).^2+PopObj(:,2).^2+PopObj(:,3).^2;
            PopCon(:,1) = (gx-9).*(4-gx);
            PopCon(:,2) = (gx-3.61).*(3.24-gx);
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,3);
            R = 1.7057*R./repmat(sqrt(sum(R.^2,2)),1,3);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            a = linspace(0,pi/2,10)';
            R = {sin(a)*cos(a')*1.7057,sin(a)*sin(a')*1.7057,cos(a)*ones(size(a'))*1.7057};
        end
    end
end
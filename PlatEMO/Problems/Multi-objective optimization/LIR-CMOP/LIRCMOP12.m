classdef LIRCMOP12 < PROBLEM
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
            obj.M = 2;
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
            sum2 = zeros(popsize,1);
            for j=2:variable_length
                if mod(j,2)==1
                    sum1=sum1+(X(:,j)-sin((0.5*j/variable_length*pi)*X(:,1))).^2;
                else
                    sum2=sum2+(X(:,j)-cos((0.5*j/variable_length*pi)*X(:,1))).^2;
                end
            end
            PopObj(:,1) = 1.7057*X(:,1).*(10*sum1+1);
            PopObj(:,2) = 1.7057*(1-X(:,1).^2).*(10*sum2+1);
            Population  = SOLUTION(X,PopObj,Constraint(PopObj),varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
           R = [1.6794,0.4419;1.3258,0.7955;0.9723,1.1490;2.0320,0.0990;
                0.6187,1.5026;0.2652,1.8562;0,2.2580;2.5690,0];
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0,5,400));
            fes   = all(Constraint([x(:),y(:)])<=0,2);  
            z     = nan(size(x));
            z(reshape(fes,size(z)) & (x/1.7057).^2+y/1.7057>=1) = 0;
            R = {x,y,z};
        end
    end
end

function PopCon = Constraint(PopObj)
    p     = 1.6;
    q     = 1.6;
    a     = 1.5;
    b     = 6;
    r     = 0.1;
    theta = -0.25 * pi;
    alpha = 0.25 * pi;
    PopCon(:,1) = r - (((PopObj(:,1)-p)*cos(theta)-(PopObj(:,2)-q)*sin(theta)).^2)/(a^2)-...
                  (((PopObj(:,1)-p)*sin(theta)+(PopObj(:,2)-q)*cos(theta)).^2)/(b^2);
    PopCon(:,2) = 2.5 - PopObj(:,1)*sin(alpha) - PopObj(:,2)*cos(alpha) + sin(4*pi*(PopObj(:,1)*cos(alpha)-PopObj(:,2)*sin(alpha)));
end
classdef LIRCMOP8 < PROBLEM
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
            for j = 2 : variable_length
                if mod(j,2) == 1
                    sum1 = sum1+(X(:,j)-sin((0.5*j/variable_length*pi)*X(:,1))).^2;
                else
                    sum2 = sum2+(X(:,j)-cos((0.5*j/variable_length*pi)*X(:,1))).^2;
                end
            end
            gx          = 0.7057;
            PopObj(:,1) = X(:,1)+10*sum1+gx;
            PopObj(:,2) = 1-X(:,1).^2+10.*sum2+gx;
            Population  = SOLUTION(X,PopObj,Constraint(PopObj),varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - sqrt(R(:,1));
            R      = R + 0.7057;
            theta  = -0.25*pi;
            c1     = 0.1 - ((R(:,1)-1.2)*cos(theta)-(R(:,2)-1.2)*sin(theta)).^2/(2^2) -...
                     ((R(:,1)-1.2)*sin(theta)+(R(:,2)-1.2)*cos(theta)).^2/(6^2);
            invalid = c1>0;
            while any(invalid)
                R(invalid,:) = (R(invalid,:)-0.7057).*1.001 + 0.7057;
                c1 = 0.1 - ((R(:,1)-1.2)*cos(theta)-(R(:,2)-1.2)*sin(theta)).^2/(2^2) -...
                     ((R(:,1)-1.2)*sin(theta)+(R(:,2)-1.2)*cos(theta)).^2/(6^2);
                invalid = c1>0;
            end
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0.7057,5,400));
            fes   = all(Constraint([x(:),y(:)])<=0,2);  
            z     = nan(size(x));
            z(reshape(fes,size(z)) & sqrt(x-0.7057)+y>=1.7057) = 0;
            R = {x,y,z};
        end
    end
end

function PopCon = Constraint(PopObj)
    p     = [1.2,2.25,3.5];
    q     = [1.2,2.25,3.5];
    a     = [2,2.5,2.5];
    b     = [6,12,10];
    r     = 0.1;
    theta = -0.25*pi;
    for k = 1 : 3
        PopCon(:,k) = r - ((PopObj(:,1)-p(k))*cos(theta)-(PopObj(:,2)-q(k))*sin(theta)).^2/(a(k)^2) -...
                      ((PopObj(:,1)-p(k))*sin(theta)+(PopObj(:,2)-q(k))*cos(theta)).^2/(b(k)^2);
    end
end
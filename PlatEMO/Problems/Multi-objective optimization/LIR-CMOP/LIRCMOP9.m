classdef LIRCMOP9 < PROBLEM
% <multi> <real> <large/none> <constrained>
% Constrained benchmark MOP with large infeasible regions

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Huang, Y. Fang, Y. You, J. Mo, C. Wei, and E.
% Goodman, An improved epsilon constraint-handling method in MOEA/D for
% CMOPs with large infeasible regions, Soft Computing, 2019, 23:
% 12491-12510.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
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
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            variable_length = size(X,2);
            popsize         = size(X,1);
            sum1            = zeros(popsize,1);
            sum2            = zeros(popsize,1);
            for j = 2 : variable_length
                if mod(j,2) == 1
                    sum1 = sum1+(X(:,j)-sin((0.5*j/variable_length*pi)*X(:,1))).^2;
                else
                    sum2 = sum2+(X(:,j)-cos((0.5*j/variable_length*pi)*X(:,1))).^2;
                end
            end
            PopObj(:,1) = 1.7057*X(:,1).*(10*sum1+1);
            PopObj(:,2) = 1.7057*(1-X(:,1).^2).*(10*sum2+1);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopCon = Constraint(obj.CalObj(X));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - R(:,1).^2;
            R      = R*1.7057;
            R(any(Constraint(R)>0,2),:) = [];
            R = [R;0,2.182;1.856,0];
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
    p     = 1.4;
    q     = 1.4;
    a     = 1.5;
    b     = 6;
    r     = 0.1;
    theta = -0.25 * pi;
    alpha = 0.25 * pi;
    PopCon(:,1) = r - (((PopObj(:,1)-p)*cos(theta)-(PopObj(:,2)-q)*sin(theta)).^2)/(a^2)-...
                  (((PopObj(:,1)-p)*sin(theta)+(PopObj(:,2)-q)*cos(theta)).^2)/(b^2);
    PopCon(:,2) = 2 - PopObj(:,1)*sin(alpha) - PopObj(:,2)*cos(alpha) + sin(4*pi*(PopObj(:,1)*cos(alpha)-PopObj(:,2)*sin(alpha)));
end
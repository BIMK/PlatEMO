classdef DASCMOP3 < PROBLEM
% <multi> <real> <large/none> <constrained>
% Difficulty-adjustable and scalable constrained benchmark MOP

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Li, C. Wei, Q. Zhang, K. Deb, and E. Goodman,
% Difficulty adjustable and scalable constrained multi-objective test
% problem toolkit, Evolutionary Computation, 2020, 28(3): 339-378.
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
        %% Calculate objective values
        function PopObj = CalObj(obj,X)       
            x_all = X(:,1:1:end); 
            sum1  = sum((x_all - sin(0.5 * pi * X(:,1))).^2,2);
            PopObj(:,1) = X(:,1) + sum1;
            PopObj(:,2) = 1.0 - sqrt(X(:,1)) + 0.5 * abs(sin(5.0 * pi * X(:,1))) + sum1;
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            x_all  = X(:,1:1:end); 
            sum1   = sum((x_all - sin(0.5 * pi * X(:,1))).^2,2);
            PopCon = Constraint(X(:,1),sum1);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,~)
            R = [0.5000,1.5000;0.5010,1.4762;0.5020,1.4710;0.5030,1.4688;0.5040,1.4681;0.6502,1.4652;0.7002,1.0541;
                 0.9044,0.8986;1.1066,0.7729;1.3008,0.6114;1.5000,0.5000;0.9069,0.8951;1.1126,0.7727;0.9129,0.8950;
                 1.1151,0.7690;0.9153,0.8914;1.1175,0.7653;1.1200,0.7616;0.9213,0.8913;1.1260,0.7613;1.1285,0.7576];
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0.5,2.5,400));
            X1    = ones(size(x));
            E     = inf(size(x));
            for a = 0 : 0.0001 : 1
                E1 = abs(a-1+sqrt(a)-0.5*abs(sin(5*pi*a))-x+y);
                X1(E1<E) = a;
                E = min(E,E1);
            end
            sum1 = x - X1;
            fes  = all(Constraint(X1(:),sum1(:))<=0,2);  
            z    = nan(size(x));
            z(reshape(fes,size(z))) = 0;
            R = {x,y,z};
        end
    end
end

function PopCon = Constraint(X1,sum1)
    % set the parameters of constraints
    DifficultyFactors = [0.5,0.5,0.5];
    % Type-I parameters
    a = 20;
    b = 2 * DifficultyFactors(1) - 1;
    % Type-II parameters
    d = 0.5;
    if DifficultyFactors(2) == 0.0
        d = 0.0;
    end
    e = d - log(DifficultyFactors(2));
    if isfinite(e) == 0
        e = 1e+30;
    end
    % Type-III parameters
    r = 0.5 * DifficultyFactors(3);
    % Calculate objective values
    PopObj(:,1) = X1 + sum1;
	PopObj(:,2) = 1.0 - sqrt(X1) + 0.5 * abs(sin(5.0 * pi * X1)) + sum1;
    % Type-I constraints
    PopCon(:,1) = b - sin(a * pi * X1);
    % Type-II constraints                               
    PopCon(:,2) = -(e - sum1) .* (sum1 - d);
    if DifficultyFactors(2) == 1.0
        PopCon(:,2) = 1e-4 - abs(sum1 - e);
    end
    % Type-III constraints
    p_k = [0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 3.0];
    q_k = [1.5, 0.5, 2.5, 1.5, 0.5, 3.5, 2.5, 1.5, 0.5];
    a_k = 0.3;
    b_k = 1.2;
    theta_k = -0.25 * pi;
    for k=1:length(p_k)
        PopCon(:,2+k) = r - ((PopObj(:,1) - p_k(k)) * cos(theta_k) - (PopObj(:,2) - q_k(k)) * sin(theta_k)).^2 ./ a_k -...
            ((PopObj(:,1) - p_k(k)) * sin(theta_k) + (PopObj(:,2) - q_k(k)) * cos(theta_k)).^2 ./ b_k;
    end
end
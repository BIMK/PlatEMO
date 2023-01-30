classdef RWMOP2 < PROBLEM
% <multi> <real> <constrained>
% Vibrating platform

%------------------------------- Reference --------------------------------
% A. Kumar, G. Wu, M. Ali, Q. Luo, R. Mallipeddi, P. Suganthan, and S. Das,
% A benchmark-suite of real-world constrained multi-objective optimization
% problems and some baseline results, Swarm and Evolutionary Computation,
% 2021, 67: 100961.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function Setting(obj)
            obj.M        = 2;
            obj.D        = 5;
            obj.lower    = [0.05,0.2,0.2,0.35,3];
            obj.upper    = [0.5,0.5,0.6,0.5,6];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            d1 = x(:,1);
            d2 = x(:,2);
            d3 = x(:,3);
            b  = x(:,4);
            L  = x(:,5);
            rho1 = 100;rho2 = 2770; rho3 = 7780;
            E1 = 1.6; E2 = 70; E3 = 200;
            c1 = 500; c2 = 1500; c3 = 800;
            mu = 2*b.*(rho1.*d1+rho2.*(d2-d1)+rho3.*(d3-d2));
            EI = (2*b./3).*(E1.*d1.^3+E2.*(d2.^3-d1.^3)+rho3.*(d3-d2));
            % Objective function
            f(:,1) = (-pi)./(2*L).^2.*(abs(EI./mu)).^0.5;
            f(:,2) = 2*b.*L.*(c1.*d1+c2.*(d2-d1)+c3.*(d3-d2));
            % Constraints
            g(:,1) = mu.*L -2800;
            g(:,2) = d1-d2;
            g(:,3) = d2-d1-0.15;
            g(:,4) = d2-d3;
            g(:,5) = d3-d2-0.01;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end 
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [-1.2746083e-03   3.1825489e+02];
        end
    end
end
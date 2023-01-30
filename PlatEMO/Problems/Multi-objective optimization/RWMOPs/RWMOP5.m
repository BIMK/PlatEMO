classdef RWMOP5 < PROBLEM
% <multi> <real> <constrained>
% Disc brake design problem

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
            obj.D        = 4;
            obj.lower    = [55,75,1000,11];
            obj.upper    = [80,110,3000,20];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            x3 = x(:,3);
            x4 = x(:,4);
            % Objective function 
            f(:,1) = 4.9e-5.*(x2.^2-x1.^2).*(x4-1);
            f(:,2) = 9.82e6.*(x2.^2-x1.^2)./(x3.*x4.*(x2.^3-x1.^3));
            % Constraints
            g(:,1) = 20 - (x2 - x1);
            g(:,2) = x3./(3.14.*(x2.^2-x1.^2)) - 0.4;
            g(:,3) = 2.22e-3.*x3.*(x2.^3-x1.^3)./(x2.^2-x1.^2).^2-1;
            g(:,4) = 900 - 2.66e-2.*x3.*x4.*(x2.^3-x1.^3)./(x2.^2-x1.^2);
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [5.3067000e+00   3.0281682e+00];
        end
    end
end
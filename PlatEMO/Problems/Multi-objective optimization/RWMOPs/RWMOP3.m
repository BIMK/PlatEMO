classdef RWMOP3 < PROBLEM
% <multi> <real> <constrained>
% Two bar truss design problem

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
            obj.D        = 3;
            obj.lower    = [1e-5,1e-5,1];
            obj.upper    = [100, 100, 3];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            x3 = x(:,3);
            % Objective function
            f(:,1) = x1.*(16+x3.^2).^(0.5)+x2.*(1+x3.^2).^(0.5);
            f(:,2) = (20.*(16+x3.^2).^(0.5))./(x3.*x1);
            % Constraints
            g(:,1) = f(:,1)-0.1;
            g(:,2) = f(:,2)-1e5;
            g(:,3) = (80.*(1+x3.^2).^(0.5))./(x3.*x2)-1e5;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [1.0000000e-01   1.0000000e+05];
        end
    end
end
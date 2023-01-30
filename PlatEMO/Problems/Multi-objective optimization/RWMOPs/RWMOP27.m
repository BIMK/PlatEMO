classdef RWMOP27 < PROBLEM
% <multi> <real> <constrained>
% Process flow sheeting problem

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
            obj.lower    = [0.2,-2.22554,-0.49];
            obj.upper    = [1,-1,1.49];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1); x2 = x(:,2); x3 = x(:,3);
            % Objective function
            f(:,1) = -0.7.*x3 + 0.8 + 5.*(0.5 - x1).^2;
            f(:,2) = x1 - x3;
            % Constraints
            g(:,1) = -(exp(x1 - 0.2) + x2);
            g(:,2) = x2 + 1.1.*x3 - 1;
            g(:,3) = x1 - x3 - 0.2;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [-2.4300000e-01   0.0000000e+00];
        end
    end
end
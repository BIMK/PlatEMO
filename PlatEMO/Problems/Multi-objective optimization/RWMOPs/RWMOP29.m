classdef RWMOP29 < PROBLEM
% <multi> <real> <constrained>
% Process synthesis problem

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
            obj.D        = 7;
            obj.lower    = [0,0,0,-0.49,-0.49,-0.49,-0.49];
            obj.upper    = [100,100,100,1.49,1.49,1.49,1.49];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1); x2 = x(:,2); x3 = x(:,3); x4 = round(x(:,4));
            x5 = round(x(:,5)); x6 = round(x(:,6)); x7 = round(x(:,7));
            % Objective function
            f(:,1) = (1-x4).^2 + (1-x5).^2 + (1-x6).^2 - log(abs(1+x7)+1e-6);
            f(:,2) = (1-x1).^2 + (2-x2).^2 + (3-x3).^2;
            % Constraints
            g(:,1) = x1 + x2 + x3 + x4 + x5 + x6 -5;
            g(:,2) = x6.^3 + x1.^2 + x2.^2 + x3.^2 - 5.5;
            g(:,3) = x1 + x4 - 1.2;
            g(:,4) = x2 + x5 - 1.8;
            g(:,5) = x3 + x6 - 2.5;
            g(:,6) = x1 + x7 - 1.2;
            g(:,7) = x5.^2 + x2.^2 - 1.64;
            g(:,8) = x6.^2 + x3.^2 - 4.25;
            g(:,9) = x5.^2 + x3.^2 - 4.64;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [2.9999990e+00   1.4000000e+01];
        end
    end
end
classdef RWMOP16 < PROBLEM
% <multi> <real> <constrained>
% Cantilever beam design problem

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
            obj.D        = 2;
            obj.lower    = [0.01 0.20];
            obj.upper    = [0.05 1];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            P  = 1;
            E  = 207000000;
            Sy = 300000;
            delta_max = 0.005;
            rho = 7800;
            % Objectives
            f(:,1) = 0.25 .* rho .* pi .* x2 .* x1.^2;
            f(:,2) = (64 .* P .* x2.^3)./(3 .* E .* pi .* x1.^4);
            % Constraints
            g(:,1) = -Sy + (32 .* P .* x2)./(pi .* x1.^3);
            g(:,2) = -delta_max + (64 .* P .* x2.^3)./(3 .* E .* pi .* x1.^4);
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [3.0630528e+00   2.0408763e-03];
        end
    end
end
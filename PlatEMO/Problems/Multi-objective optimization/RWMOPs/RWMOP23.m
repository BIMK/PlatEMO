classdef RWMOP23 < PROBLEM
% <multi> <real> <constrained>
% Reactor network design

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
            obj.D        = 6;
            obj.lower    = [0,0,0,0,0.00001,0.00001];
            obj.upper    = [1,1,1,1,16,16];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            k1 = 0.09755988; k2 = 0.99*k1;
            k3 = 0.0391908; k4 = 0.9*k3;
            x1 = x(:,1); x2 = x(:,2); x3 = x(:,3);
            x4 = x(:,4); x5 = x(:,5); x6 = x(:,6);
            % Objective function
            f(:,1) = -x4;
            f(:,2) = x5.^(0.5)+x6.^(0.5);
            % Constraints
            g(:,1) = f(:,2)-4;
            h(:,1) = k1.*x5.*x2 + x1 -1;
            h(:,2) = k3.*x5.*x3+x3+x1-1;
            h(:,3) = k2.*x6.*x2 - x1 + x2;
            h(:,4) = k4.*x6.*x4 + x2-x1+x4-x3;
            g      = [g,abs(h)-1e-4];
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [-4.0194083e-04   4.0000000e+00];
        end
    end
end
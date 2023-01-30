classdef RWMOP28 < PROBLEM
% <multi> <real> <constrained>
% Two reactor problem

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
            obj.lower    = [0,0,0,0,-0.49,-0.49,0];
            obj.upper    = [20,20,10,10,1.49,1.49,40];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            v1 = x(:,3);
            v2 = x(:,4);
            y1 = round(x(:,5));
            y2 = round(x(:,6));
            x_ = x(:,7);
            z1 = 0.9*(1-exp(-0.5.*v1)).*x1;
            z2 = 0.8*(1-exp(-0.4*v2)).*x2;
            % Objective function
            f(:,1) = 7.5.*y1 + 5.5.*y2 + 7.*v1 + 6.*v2 + 5.*x_;
            f(:,2) = x1 + x2;
            % Constraints
            h(:,1) = y1 + y2 - 1;
            h(:,2) = z1 + z2 - 10;
            h(:,3) = x1 + x2 -x_;
            h(:,4) = z1.*y1 + z2.*y2 - 10;
            g(:,1) = v1 - 10*y1 -1e-6;
            g(:,2) = v2 - 10*y2;
            g(:,3) = x1 - 20*y1;
            g(:,4) = x2 - 20*y2;
            g      = [g,abs(h)-1e-4];
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [1.3343214e+02   1.3108200e+01];
        end
    end
end
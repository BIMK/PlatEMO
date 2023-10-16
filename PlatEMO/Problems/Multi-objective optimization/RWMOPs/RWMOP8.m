classdef RWMOP8 < PROBLEM
% <multi> <real> <constrained>
% Car side impact design problem

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
            obj.M        = 3;
            obj.D        = 7;
            obj.lower    = [0.5,0.45,0.5,0.5,0.875,0.4,0.4];
            obj.upper    = [1.5,1.35,1.5,1.5,2.625,1.2,1.2];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            x3 = x(:,3);
            x4 = x(:,4);
            x5 = x(:,5);
            x6 = x(:,6);
            x7 = x(:,7);
            VMBP = 10.58-0.674.*x1.*x2-0.67275.*x2;
            VFD  = 16.45-0.489.*x3.*x7-0.843.*x5.*x6;
            % Objective function
            f(:,1) = 1.98+4.9.*x1.*6.67.*x2+6.98.*x3+4.01.*x4+1.78.*x5+1e-5.*x6+2.73.*x7;
            f(:,2) = 4.72-0.5.*x4-0.19.*x2.*x3;
            f(:,3) = 0.5.*(VMBP+VFD);
            % Constraints
            g(:,1) = -1+1.16-0.3717.*x2.*x4-0.0092928.*x3;
            g(:,2) = -0.32+0.261-0.0159.*x1.*x2-0.06486.*x1-0.019.*x2.*x7+0.0144.*x2.*x5+0.0154464.*x6;
            g(:,3) = -0.32+0.74-0.61.*x2-0.031296.*x3-0.031872.*x7+0.227.*x2.^2;
            g(:,4) = -0.32+0.214+0.00817.*x5-0.045195.*x1-0.0135168.*x1+0.03099.*x2.*x6-0.018.*x2.*x7+0.007176.*x3+0.023232.*x3-0.00364.*x5.*x6-0.018.*x2.^2;
            g(:,5) = -32+33.86+2.95.*x3-5.057.*x1.*x2-3.795.*x2-3.4431.*x7+1.45728;
            g(:,6) = -32+28.98+3.818.*x3-4.2.*x1.*x2+1.27296.*x6-2.68065.*x7;
            g(:,7) = -32+46.36-9.9.*x2-4.4505.*x1;
            g(:,8) = f(:,2)-4;
            g(:,9) = VMBP - 9.9;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [9.2596587e+01   4.0000000e+00   1.2699733e+01];
        end
    end
end
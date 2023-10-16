classdef RWMOP7 < PROBLEM
% <multi> <real> <constrained>
% Gear train design problem

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
            obj.lower    = [11.51,11.51,11.51,11.51];
            obj.upper    = [60.49,60.49,60.49,60.49];
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
            f(:,1) = abs(6.931-x3.*x4./(x1.*x2));
            f(:,2) = max(x,[],2);
            % Constraints
            g(:,1) = f(:,1)./6.931-0.5;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [3.4655000e+00   4.5417795e+01];
        end
    end
end
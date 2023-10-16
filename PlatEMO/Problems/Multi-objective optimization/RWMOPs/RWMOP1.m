classdef RWMOP1 < PROBLEM
% <multi> <real> <constrained>
% Pressure vessal problem

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
            obj.lower    = [0.51,0.51,10,10];
            obj.upper    = [99.49,99.49,200,200];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = round(x(:,1));
            x2 = round(x(:,2));
            x3 = x(:,3);
            x4 = x(:,4);
            z1 = 0.0625*x1;
            z2 = 0.0625*x2;
            % Objective function
            f(:,1) = 1.7781.*z1.*x3.^2+0.6224.*z1.*x2.*x4+3.1661.*z1.^2.*x4+19.84.*z1.^2.*x3;
            f(:,2) = -pi.*x3.^2.*x4-(4/3).*pi.*x3.^3;
            % Constraints
            g(:,1) = 0.00954.*x3-z2;
            g(:,2) = 0.0193.*x3-z1;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [3.5964885e+05  -7.3303829e+03];
        end
    end
end
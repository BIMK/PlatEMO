classdef RWMOP22 < PROBLEM
% <multi> <real> <constrained>
% Haverly's pooling problem

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
            obj.D        = 9;
            obj.lower    = zeros(1,9);
            obj.upper    = [100,200,100,100,100,100,200,100,200];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1); x2 = x(:,2); x3 = x(:,3);
            x4 = x(:,4); x5 = x(:,5); x6 = x(:,6);
            x7 = x(:,7); x8 = x(:,8); x9 = x(:,9);
            % Objective function
            f(:,1) = -9*x1-15*x2+6*x3+16*x4;
            f(:,2) = 10.*(x5+x6);
            % Constraints
            g(:,1) = x9.*x7+2*x5-2.5*x1;
            g(:,2) = x9.*x8+2*x6-1.5*x2;
            h(:,1) = x7+x8-x4-x3;
            h(:,2) = x1-x5-x7;
            h(:,3) = x2-x6-x8;
            h(:,4) = x9.*x7+x9.*x8-3.*x3-x4;
            g      = [g,abs(h)-1e-4];
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [-1.0590871e+02   2.0000000e+03];
        end
    end
end
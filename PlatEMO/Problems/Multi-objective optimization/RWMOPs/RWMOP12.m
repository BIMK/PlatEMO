classdef RWMOP12 < PROBLEM
% <multi> <real> <constrained>
% Simply supported I-beam design

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
            obj.lower    = [10 10 0.9 0.9];
            obj.upper    = [80 50 5 5];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            x3 = x(:,3);
            x4 = x(:,4);
            P  = 600;
            L  = 200;
            E  = 2e4;
            % Objectives
            f(:,1) = 2 .* x2 .* x4 + x3 .* (x1- 2.*x4);
            f(:,2) = P.*L.^3./(48.*E.*(x3 .*((x1 - 2.*x4).^3)+2.*x2.*x4.*(4.*x4.*x4+3.*x1.*(x1-2.*x4)))./12);
            % Constraints
            g(:,1) = -16 + 180000*x1./(x3.*((x1 - 2*x4).^3)+2*x2.*x4.*(4.*x4.*x4+3*x1.*(x1-2*x4))) + 15000.*x2./((x1-2.*x4).*x3.^3+2.*x4.*x2.^3);
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [4.3795472e+02   6.1459097e-02];
        end
    end
end
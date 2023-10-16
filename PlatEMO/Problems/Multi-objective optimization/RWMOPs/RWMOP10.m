classdef RWMOP10 < PROBLEM
% <multi> <real> <constrained>
% Two bar plane truss

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
            obj.lower    = [0.1,0.5];
            obj.upper    = [2,2.5];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x    = varargin{1};
            x1   = x(:,1);
            x2   = x(:,2);
            rho  = 0.283;
            h    = 100;
            P    = 104;
            E    = 3e7;
            rho0 = 2e4;
            % Objective function
            f(:,1) = 2.*rho.*h.*x2.*sqrt(1+x1.^2);
            f(:,2) = rho.*h.*(1+x1.^2).^1.5.*(1+x1.^4).^0.5./(2.*sqrt(2).*E.*x1.^2.*x2);
            % Constraints
            g(:,1) = P.*(1+x1).*(1+x1.^2).^0.5./(2.*sqrt(2).*x1.*x2)-rho0;
            g(:,2) = P.*(-x1+1).*(1+x1.^2).^0.5./(2.*sqrt(2).*x1.*x2)-rho0;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [1.8704862e+02   6.7710178e-05];
        end
    end
end
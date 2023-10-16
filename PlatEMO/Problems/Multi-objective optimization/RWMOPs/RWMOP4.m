classdef RWMOP4 < PROBLEM
% <multi> <real> <constrained>
% Weldan beam design problem

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
            obj.lower    = [0.125,0.1,0.1,0.125];
            obj.upper    = [5,10,10,5];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            x3 = x(:,3);
            x4 = x(:,4);
            P  = 6000;
            L  = 14;
            E  = 30e6;
            tmax   = 13600;
            sigmax = 30000;
            G      = 12e6;
            Pc     = (4.013.*E.*((x3.^2+x4.^6)./36).^0.5)./(L.^2).*(1-x3./(2.*L).*(E./(4*G))^(0.5));
            sigma  = (6.*P.*L)./(x4.*x3.^2);
            J      = 2.*(sqrt(2).*x1.*x2.*(x2.^2/12+((x1+x3)./2).^2));
            R      = sqrt(x2.^2/4+((x1+x3)./2).^2);
            M      = P.*(L+x2./2);
            tho1   = P./(sqrt(2).*x1.*x2);
            tho2   = M.*R./J;
            tho    = sqrt(tho1.^2+2.*tho1.*tho2.*x2./(2*R)+tho2.^2);
            % Objective function
            f(:,1) = 1.10471.*x1.^2.*x2+0.04811.*x3.*x4.*(14+x2);
            f(:,2) = (4.*P.*L.^3)./(E.*x4.*x3.^3);
            % Constraints
            g(:,1) = tho - tmax;
            g(:,2) = sigma - sigmax;
            g(:,3) = x1 - x4;
            g(:,4) = P - Pc;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [3.6679325e+01   1.3066667e-02];
        end
    end
end
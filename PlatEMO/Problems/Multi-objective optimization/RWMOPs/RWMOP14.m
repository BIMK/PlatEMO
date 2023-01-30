classdef RWMOP14 < PROBLEM
% <multi> <real> <constrained>
% Multiple-disk clutch brake design problem

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
            obj.D        = 5;
            obj.lower    = [60 90 1 0 2];
            obj.upper    = [80 110 3 1000 9];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            Mf = 3; Ms = 40; Iz = 55; n = 250; Tmax = 15; s = 1.5; delta = 0.5; 
            Vsrmax = 10; rho = 0.0000078; pmax = 1; mu = 0.6; Lmax = 30; delR = 20;
            Rsr = 2./3.*(x(:,2).^3-x(:,1).^3)./(x(:,2).^2.*x(:,1).^2);
            Vsr = pi.*Rsr.*n./30;
            A   = pi.*(x(:,2).^2-x(:,1).^2);
            Prz = x(:,4)./A;
            w   = pi.*n./30;
            Mh  = 2/3.*mu.*x(:,4).*x(:,5).*(x(:,2).^3-x(:,1).^3)./(x(:,2).^2-x(:,1).^2);
            T   = Iz.*w./(Mh+Mf);
            % Objective function
            f(:,1) = pi.*(x(:,2).^2-x(:,1).^2).*x(:,3).*(x(:,5)+1).*rho;
            f(:,2) = T;
            % Constraints
            g(:,1) = -x(:,2)+x(:,1)+delR;
            g(:,2) = (x(:,5)+1).*(x(:,3)+delta)-Lmax;
            g(:,3) = Prz-pmax;
            g(:,4) = Prz.*Vsr-pmax.*Vsrmax;
            g(:,5) = Vsr-Vsrmax;
            g(:,6) = T-Tmax;
            g(:,7) = s.*Ms-Mh;
            g(:,8) = -T;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [1.3967521e+00   1.4920745e-02];
        end
    end
end
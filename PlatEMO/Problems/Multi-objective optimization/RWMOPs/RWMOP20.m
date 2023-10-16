classdef RWMOP20 < PROBLEM
% <multi> <real> <constrained>
% Hydro-static thrust bearing design problem

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
            obj.lower    = [ 1, 1,  1e-6,1];
            obj.upper    = [16, 16, 16*1e-6,16];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x = varargin{1};
            R = x(:,1); Ro = x(:,2);  mu = x(:,3); Q = x(:,4);
            gamma = 0.0307; C = 0.5; n = -3.55; C1 = 10.04;
            Ws = 101000; Pmax = 1000; delTmax = 50; hmin = 0.001;
            gg = 386.4; N = 750;
            P    = (log10(log10(8.122*1e6.*mu+0.8))-C1)./n;
            delT = 2.*(10.^P-560);
            Ef   = 9336.*Q.*gamma.*C.*delT;
            h    = (2.*pi.*N./60).^2.*2.*pi.*mu./Ef.*(R.^4./4-Ro.^4./4)-1e-5;
            Po   = (6.*mu.*Q./(pi.*h.^3)).*log(R./Ro);
            W    = pi.*Po./2.*(R.^2-Ro.^2)./(log(R./Ro)-1e-5);
            % Objective function
            f(:,1) = (Q.*Po./0.7+Ef)./12;
            f(:,2) = gamma./(gg.*Po).*(Q./(2.*pi.*R.*h));
            % Constraints
            g(:,1) = Ws-W;
            g(:,2) = Po-Pmax;
            g(:,3) = delT-delTmax;
            g(:,4) = hmin-h;
            g(:,5) = Ro-R;
            g(:,6) = f(:,2)-0.001;
            g(:,7) = W./(pi.*(R.^2-Ro.^2)+1e-5)-5000;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [2.6725846e+02  -2.7672651e-05];
        end
    end
end
classdef RWMOP17 < PROBLEM
% <multi> <real> <constrained>
% Bulk carriers design problem 

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
            obj.D        = 6;
            obj.lower    = [150.0 20.0 13.0 10.0 14.0 0.63];
            obj.upper    = [274.32 32.31 25.0 11.71 18.0 0.75];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x   = varargin{1};
            L   = x(:,1);
            B   = x(:,2);
            D   = x(:,3);
            T   = x(:,4);
            V_k = x(:,5);
            C_B = x(:,6);

            a   = 4977.06.*C_B.^2 - 8105.61.*C_B + 4456.51;
            b   = -10847.2.*C_B.^2 + 12817.*C_B - 6960.32;
            F_n = 0.5144./(9.8065 .* L).^0.5;
            P   = ((1.025.*L.*B.*T.*C_B).^(2/3).*V_k.^3)./(a + b.*F_n);

            W_s = 0.034.*L.^1.7.*B.^0.6.*D.^0.4.*C_B.^0.5;
            W_o = L.^0.8.*B.^0.6.*D.^0.3.*C_B.^0.1;
            W_m = 0.17.*P.^0.9;
            ls  = W_s+W_o+W_m;

            D_wt  = 1.025.*L.*B.*T.*C_B-ls;
            F_c   = 0.19.*24.*P./1000 + 0.2;
            D_cwt = D_wt - F_c.*((5000.*V_k)./24 + 5)-2.*D_wt.^0.5;
            R_trp = 350./((5000.*V_k)./24 + 2.*(D_cwt./8000 + 0.5));
            ac    = D_cwt.*R_trp;
            S_d   = 5000.*V_k./24;

            C_c = 0.2.*1.3 .* (2000.*W_s.^0.85 + 3500.*W_o + 2400.*P.^0.8);
            C_r = 40000.*D_wt.^0.3;
            C_v = (1.05.*100.*F_c.*S_d + 6.3.*D_wt.^0.8).*R_trp;
            
            % Objectives
            f(:,1) = (C_c + C_r + C_v)./ac;
            f(:,2) = ls;
            f(:,3) = -ac;
            
            % Constraints
            g(:,1) = L./B - 6;
            g(:,2) = 15 - L./D;
            g(:,3) = 19 - L./T;
            g(:,4) = 0.45.*D_wt.^0.31 - T;
            g(:,5) = 0.7.*D + 0.7 - T;
            g(:,6) = 0.32 - F_n;
            g(:,7) = 0.53.*T + ((0.085.*C_B - 0.002).*B.^2)./(T.*C_B)-(1 + 0.52.*D) - 0.07.*B;
            g(:,8) = D_wt - 3000;
            g(:,9) = 500000 - D_wt;
            g      = -g;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [-3.1514157e+03   8.2606298e+03   8.1260004e+02];
        end
    end
end
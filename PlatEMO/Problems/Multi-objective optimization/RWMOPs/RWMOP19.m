classdef RWMOP19 < PROBLEM
% <multi> <real> <constrained>
% Multi-product batch plant

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
            obj.D        = 10;
            obj.lower    = [0.51,0.51,0.51,250,250,250,6,4,40,10];
            obj.upper    = [3.49,3.49,3.49,2500,2500,2500,20,16,700,450];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x = varargin{1};
            S = [2,3,4;
                4,6,3];
            t = [8,20,8;
                16,4,4];
            H = 6000; alp = 250; beta = 0.6;
            Q1 = 40000; Q2 = 20000;
            N1 = round(x(:,1)); N2 = round(x(:,2)); N3 = round(x(:,3));
            V1 = x(:,4); V2 = x(:,5); V3 = x(:,6);
            TL1 = x(:,7); TL2 = x(:,8);
            B1 = x(:,9); B2 = x(:,10);
            % Objective function
            f(:,1) = alp.*(N1.*V1.^beta+N2.*V2.^beta+N3.*V3.^beta);
            f(:,2) = 65.*(Q1./B1+Q2./B2)+0.08.*Q1+0.1.*Q2;
            f(:,3) = Q1.*TL1./B1+Q2.*TL2./B2;
            % Constraints
            g(:,1) = Q1.*TL1./B1+Q2.*TL2./B2-H;
            g(:,2) = S(1,1).*B1+S(2,1).*B2-V1;
            g(:,3) = S(1,2).*B1+S(2,2).*B2-V2;
            g(:,4) = S(1,3).*B1+S(2,3).*B2-V3;
            g(:,5) = t(1,1)-N1.*TL1;
            g(:,6) = t(1,2)-N2.*TL1;
            g(:,7) = t(1,3)-N3.*TL1;
            g(:,8) = t(2,1)-N1.*TL2;
            g(:,9) = t(2,2)-N2.*TL2;
            g(:,10) = t(2,3)-N3.*TL2;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [2.4435182e+05  4.8572551e+04   6.0000000e+03];
        end
    end
end
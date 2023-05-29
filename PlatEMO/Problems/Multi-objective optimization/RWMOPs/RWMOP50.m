classdef RWMOP50 < PROBLEM
% <multi> <real> <constrained>
% Power distribution system planning

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
            obj.D        = 6;
            obj.lower    = [10 10 35 35 125 130];
            obj.upper    = [125 150 210 225 315 325];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x = varargin{1};
            PD = 1200;
            B = [140 17 15 19 26 22;
                  17 60 13 16 15 20;
                  15 13 65 17 24 19;
                  19 16 17 71 30 25;
                  26 15 24 30 69 32;
                  22 20 19 25 32 85] * 10^-6;
            a     = [756.7988 451.3251 1243.5311 1049.9977 1356.6592 1658.5696];
            b     = [38.5390  46.1591 38.3055  40.3965 38.2704 36.3278];
            c     = [0.15247  0.10587  0.03546  0.02803  0.01799 0.02111];
            alpha = [13.8593 13.8593  40.2669  40.2669 42.8955 42.8955];
            beta  = [0.32767  0.32767 -0.54551 -0.54551 -0.51116 -0.51116];
            gamma = [0.00419 0.00419  0.00683  0.00683 0.00461 0.00461];
            PL = zeros(size(x,1),1);
            for i = 1 : size(x,2)
                for j = 1 : size(x,2)
                    PL = PL + x(:,i) .* B(i,j) .* x(:,j);
                end
            end
            % Objectives
            f = zeros(size(x,1),2);
            for i = 1 : size(x,2)
                f(:,1) = f(:,1) + a(i) + b(i) .* x(:,i) + c(i) .* x(:,i).^2 ;
            end
            for i = 1 : size(x,2)
                f(:,2) = f(:,2) + alpha(i) + beta(i) .* x(:,i) + gamma(i) .* x(:,i).^2;
            end
            % Constraints
            h = sum(x,2) - PD - PL;
            Population = SOLUTION(varargin{1},f,abs(h)-1e-4,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [6.4895619e+04   1.2538936e+03];
        end
    end
end
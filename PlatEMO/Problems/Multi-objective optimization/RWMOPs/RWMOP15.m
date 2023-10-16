classdef RWMOP15 < PROBLEM
% <multi> <real> <constrained>
% Spring design problem

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
            obj.D        = 3;
            obj.lower    = [0.51,0.6,0.51];
            obj.upper    = [70.49,3,42.49];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = round(x(:,1));
            x2 = x(:,2);
            d  = [0.009,0.0095,0.0104,0.0118,0.0128,0.0132,0.014,....
                  0.015, 0.0162, 0.0173, 0.018, 0.020, 0.023, 0.025,...
                  0.028, 0.032, 0.035, 0.041, 0.047, 0.054, 0.063,....
                  0.072, 0.080, 0.092, 0.0105, 0.120, 0.135, 0.148,....
                  0.162, 0.177, 0.192, 0.207, 0.225, 0.244, 0.263,....
                  0.283, 0.307, 0.331, 0.362,0.394,0.4375,0.500];
            x3 = d(max(1,min(42,round(x(:,3))))); x3 = x3(:);
            cf   = (4.*x2./x3-1)./(4.*x2./x3-4)+0.615.*x3./x2;
            K    = (11.5.*10.^6.*x3.^4)./(8.*x1.*x2.^3);
            lf   = 1000./K + 1.05.*(x1+2).*x3;
            sigp = 300./K;
            % Objective function
            f(:,1) = (pi.^2.*x2.*x3.^2.*(x1+2))./4;
            f(:,2) = (8000.*cf.*x2)./(pi.*x3.^3);
            % Constraints
            g(:,1) = (8000.*cf.*x2)./(pi.*x3.^3)-189000;
            g(:,2) = lf-14;
            g(:,3) = 0.2-x3;
            g(:,4) = x2-3;
            g(:,5) = 3-x2./x3;
            g(:,6) = sigp - 6;
            g(:,7) = sigp+700./K+1.05.*(x1+2).*x3-lf;
            g(:,8) = 1.25-700./K;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [2.7942662e+01   1.8799119e+05];
        end
    end
end
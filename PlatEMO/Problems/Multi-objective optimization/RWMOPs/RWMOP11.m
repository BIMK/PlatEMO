classdef RWMOP11 < PROBLEM
% <multi> <real> <constrained>
% Water resource management problem

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
            obj.M        = 5;
            obj.D        = 3;
            obj.lower    = [0.01 0.01 0.01];
            obj.upper    = [0.45 0.1 0.1];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            x3 = x(:,3);
            % Objectives
            f(:,1) = 106780.37 .* (x2 + x3) + 61704.67 ;
            f(:,2) = 3000 .* x1 ;
            f(:,3) = 305700 .* 2289 .* x2 ./ power(0.06.*2289, 0.65) ;
            f(:,4) = 250 .* 2289 .* exp(-39.75.*x2+9.9.*x3+2.74) ;
            f(:,5) = 25 .* (1.39 ./(x1.*x2) + 4940.*x3 -80) ;
            % Constraints   
            g(:,1) = 1 - (0.00139./(x1.*x2)+4.94.*x3-0.08);
            g(:,2) = 1 - (0.000306./(x1.*x2)+1.082.*x3-0.0986);
            g(:,3) = 50000 - (12.307./(x1.*x2) + 49408.24.*x3+4051.02);
            g(:,4) = 16000 - (2.098./(x1.*x2)+8046.33.*x3-696.71);
            g(:,5) = 10000 - (2.138./(x1.*x2)+7883.39.*x3-705.04);
            g(:,6) = 2000 - (0.417.*x1.*x2 + 1721.26.*x3-136.54);
            g(:,7) = 550 - (0.164./(x1.*x2)+631.13.*x3-54.48);
            g      = -g;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [7.3450511e+04   1.3500000e+03   2.8534690e+06   6.6200320e+06   2.5000000e+04];
        end
    end
end
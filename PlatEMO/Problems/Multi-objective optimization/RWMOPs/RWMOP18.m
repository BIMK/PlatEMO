classdef RWMOP18 < PROBLEM
% <multi> <real> <constrained>
% Front rail design problem

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
            obj.lower    = [136,56,1.4];
            obj.upper    = [146,68,2.2];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            hh = x(:,1);
            w  = x(:,2);
            t  = x(:,3);

            Ea = 14496.5;
            Fa = 234.9;
            E  = -70973.4 + 958.656 .* w+ 614.173 .* hh - 3.827 .* w .* hh + 57.023 .* w .* t + 63.274 .* hh .* t-3.582 .* w.^2 - 1.4842 .* hh.^2 - 1890.174 .* t.^2;
            F  = 111.854 - 20.210 .* w + 7.560 .* hh - 0.025 .* w .* hh + 2.731 .* w .* t - 1.479 .* hh .* t + 0.165 .* w.^2;

            % Objectives
            f(:,1) = Ea./E;
            f(:,2) = F./Fa;

            % Constraints
            g(:,1) = (hh - 136).* (146 - hh);
            g(:,2) = (w - 58) .* (66 - w);
            g(:,3) = (t - 1.4) .* (2.2 - t);
            g      = -g;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [9.3366407e-01   1.1965960e+00];
        end
    end
end
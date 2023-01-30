classdef RWMOP31 < PROBLEM
% <multi> <real> <constrained>
% Synchronous pptimal pulse-width modulation of 5-level inverters

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
            obj.D        = 25;
            obj.lower    = zeros(1,25);
            obj.upper    = 90*ones(1,25);
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x = varargin{1};
            m = 0.32;
            s = [1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,1,1,-1];
            k = [5,7,11,13,17,19,23,25,29,31,35,37,41,43,47,49,53,55,59,61,65,67,71,73,77,79,83,85,91,95,97];
            % Objective function
            for i = 1 : size(x,1)
                su = 0;
                for j = 1 : 31
                    su2 = 0;
                    for l = 1 : size(x,2)
                        su2 = su2 + s(l).*cos(k(j).*x(i,l)*pi/180);
                    end
                    su = su + su2.^2./k(j).^4;
                end
                f(i,1) = (su).^0.5./(sum(1./k.^4)).^0.5;
            end
            f(:,2) = (sum(s.*cos(x*pi/180),2)-m).^2;
            % Constraints
            for i = 1 : size(x,2)-1
                g(:,i) = x(:,i)-x(:,i+1)+1e-6;
            end
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [7.3429846e-01   1.0239997e-01];
        end
    end
end
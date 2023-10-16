classdef RWMOP6 < PROBLEM
% <multi> <real> <constrained>
% Speed reducer design problem

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
            obj.D        = 7;
            obj.lower    = [2.6,0.7,16.51,7.3,7.3,2.9,5];
            obj.upper    = [3.6,0.8,28.49,8.3,8.3,3.9,5.5];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1);
            x2 = x(:,2);
            x3 = round(x(:,3));
            x4 = x(:,4);
            x5 = x(:,5);
            x6 = x(:,6);
            x7 = x(:,7);
            % Objective function
            f(:,1) = 0.7854.*x1.*x2.^2.*(10.*x3.^2/3+14.933.*x3-43.0934) - 1.508.*x1.*(x6.^2+x7.^2)+7.477.*(x6.^3+x7.^3)+0.7854.*(x4.*x6.^2+x5.*x7.^2);
            f(:,2) = sqrt((745.*x4./(x2.*x3)).^2+1.69e7)./(0.1.*x6.^3);
            % Constraints
            g(:,1) = 1./(x1.*x2.^2.*x3)-1/27;
            g(:,2) = 1./(x1.*x2.^2.*x3.^2)-1/397.5;
            g(:,3) = x4.^3./(x2.*x3.*x6.^4)-1/1.93;
            g(:,4) = x5.^3./(x2.*x3.*x7.^4)-1/1.93;
            g(:,5) = x2.*x3-40;
            g(:,6) = x1./x2-12;
            g(:,7) = -x1./x2+5;
            g(:,8) = 1.9-x4+1.5.*x6;
            g(:,9) = 1.9-x5+1.1.*x7;
            g(:,10) = f(:,2)-1300;
            g(:,11) = sqrt((745.*x5./(x2.*x3)).^2+1.575e8)./(0.1*x7.^3)-850;
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
         %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [5.9698221e+03   1.3000000e+03];
        end
    end
end
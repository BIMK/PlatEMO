classdef RWMOP21 < PROBLEM
% <multi> <real> <constrained>
% Crash energy management for high-speed train

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
            obj.lower    = [1.3, 2.5, 1.3, 1.3, 1.3, 1.3];
            obj.upper    = [1.7, 3.5, 1.7, 1.7, 1.7, 1.7];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x  = varargin{1};
            x1 = x(:,1); x2 = x(:,2); x3 = x(:,3);
            x4 = x(:,4); x5 = x(:,5); x6 = x(:,6);
            % Objective function
            f(:,1) = 1.3667145844797-0.00904459793976106*x1-0.0016193573938033*x2-0.00758531275221425*x3-0.00440727360327102*x4-0.00572216860791644*x5-0.00936039926190721*x6+2.62510221107328*10^(-6)*(x1.^2)+4.92982681358861*10^(-7)*(x2.^2)+2.25524989067108*10^(-6)*(x3.^2)+...
            1.84605439400301*10^(-6)*(x4.^2)+2.17175358243416*10^(-6)*(x5.^2)+3.90158043948054*10^(-6)*(x6.^2)+4.55276994245781*10^(-7)*x1.*x2-6.37013576290982*10^(-7)*x1.*x3+8.26736480446359*10^(-7)*x1.*x4+5.66352809442276*10^(-8)*x1.*x5-3.20213897443278*10^(-7)*x1.*x6+...
            1.18015467772812*10^(-8)*x2.*x3+9.25820391546515*10^(-8)*x2.*x4-1.05705364119837*10^(-7)*x2.*x5-4.74797783014687*10^(-7)*x2.*x6-5.02319867013788*10^(-7)*x3.*x4+9.54284258085225*10^(-7)*x3.*x5+1.80533309229454*10^(-7)*x3.*x6-1.07938022118477*10^(-6)*x4.*x5-...
            1.81370642220182*10^(-7)*x4.*x6-2.24238851688047*10^(-7)*x5.*x6; 
            f(:,2) = -1.19896668942683+3.04107017009774*x1+1.23535701600191*x2+2.13882039381528*x3+2.33495178382303*x4+2.68632494801975*x5+3.43918953617606*x6-7.89144544980703*10^(-4)*(x1.^2)-2.06085185698215*10^(-4)*(x2.^2)-7.15269900037858*10^(-4)*(x3.^2)-7.8449237573837*10^(-4)*(x4.^2)-...
            9.31396896237177*10^(-4)*(x5.^2)-1.40826531972195*10^(-3)*(x6.^2)-1.60434988248392*10^(-4)*x1.*x2+2.0824655419411*10^(-4)*x1.*x3-3.0530659653553*10^(-4)*x1.*x4-8.10145973591615*10^(-5)*x1.*x5+6.94728759651311*10^(-5)*x1.*x6+1.18015467772812*10^(-8)*x2.*x3+...
            9.25820391546515*10^(-8)*x2.*x4-1.05705364119837*10^(-7)*x2.*x5+1.69935290196781*10^(-4)*x2.*x6+2.32421829190088*10^(-5)*x3.*x4-2.0808624041163476*10^(-4)*x3.*x5+1.75576341867273*10^(-5)*x3.*x6+2.68422081654044*10^(-4)*x4.*x5+4.39852066801981*10^(-5)*x4.*x6+...
            2.96785446021357*10^(-5)*x5.*x6;
            % Constraints
            g(:,1) = f(:,1)-5;
            g(:,2) = -f(:,1);
            g(:,3) = f(:,2) - 28;
            g(:,4) = -f(:,2);
            Population = SOLUTION(varargin{1},f,g,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [1.3157336e+00   2.6297736e+01];
        end
    end
end
classdef DASCMOP5 < PROBLEM
% <problem> <DAS-CMOP>
% Difficulty-adjustable and scalable constrained benchmark MOP

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Li, C. Wei, Q. Zhang, K. Deb, and E. Goodman,
% Difficulty adjustable and scalable constrained multi-objective test
% problem toolkit, Evolutionary Computation, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenji Li
    
    methods
        %% Initialization
        function obj = DASCMOP5()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 30;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X) 
            x_all = X(:,2:1:end); 
            g_1   = sum((x_all - 0.5) .* (x_all - 0.5) - cos(20.0 * pi .* ( x_all - 0.5)),2);
            sum1  = obj.Global.D - 1 + g_1;
            
            PopObj(:,1) = X(:,1) + sum1;
            PopObj(:,2) = 1.0 - sqrt(X(:,1)) + sum1;
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            x_all = X(:,2:1:end); 
            g_1   = sum((x_all - 0.5) .* (x_all - 0.5) - cos(20.0 * pi .* ( x_all - 0.5)),2);
            sum1  = obj.Global.D - 1 + g_1;
            
            % set the parameters of constraints
            DifficultyFactors = [0.5,0.5,0.5];
            % Type-I parameters
            a = 20;
            b = 2 * DifficultyFactors(1) - 1;
            % Type-II parameters
            d = 0.5;
            if DifficultyFactors(2) == 0.0
                d = 0.0;
            end
            e = d - log(DifficultyFactors(2));
            if isfinite(e) == 0
                e = 1e+30;
            end
            % Type-III parameters
            r = 0.5 * DifficultyFactors(3);
            
            PopObj = obj.CalObj(X);
            % Type-I constraints
            PopCon(:,1) = b - sin(a * pi * X(:,1));
            % Type-II constraints                      
            PopCon(:,2) = -(e - sum1) .* (sum1 - d);
            if DifficultyFactors(2) == 1.0
                PopCon(:,2) = 1e-4 - abs(sum1 - e);
            end
            % Type-III constraints
            p_k = [0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 3.0];
            q_k = [1.5, 0.5, 2.5, 1.5, 0.5, 3.5, 2.5, 1.5, 0.5];
            a_k = 0.3;
            b_k = 1.2;
            theta_k = -0.25 * pi;
            for k=1:length(p_k)
                PopCon(:,2+k) = r - ((PopObj(:,1) - p_k(k)) * cos(theta_k) - (PopObj(:,2) - q_k(k)) * sin(theta_k)).^2 ./ a_k -...
                    ((PopObj(:,1) - p_k(k)) * sin(theta_k) + (PopObj(:,2) - q_k(k)) * cos(theta_k)).^2 ./ b_k;
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'DASCMOP_PF.mat'),'PF');
            P = PF{5};
        end 
    end
end
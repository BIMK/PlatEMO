classdef BBOB_F19 < PROBLEM
% <2009> <single> <real> <expensive/none>
% Composite Griewank-Rosenbrock function F8F2
% xopt --- 0 --- Optimal decision variables

%------------------------------- Reference --------------------------------
% N. Hansen, S. Finck, R. Ros, and A. Auger. Real-parameter black-box
% optimization benchmarking 2009: Noiseless functions definitions. RR-6829,
% INRIA. 2009.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        xopt;   % Optimal decision variables
        R;      % Orthogonal (rotation) matrix
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.xopt = obj.ParameterSet(0);
            obj.M = 1;
            if isscalar(obj.xopt)
                if isempty(obj.D)
                    obj.D = 30;
                end
                obj.xopt = repmat(obj.xopt,1,obj.D);
            else
                obj.D = length(obj.xopt);
            end
            obj.lower    = zeros(1,obj.D) - 5;
            obj.upper    = zeros(1,obj.D) + 5;
            obj.encoding = ones(1,obj.D);
            obj.R = orth(rand(RandStream('mlfg6331_64','Seed',2),obj.D));
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            Z = max(1,sqrt(obj.D)/8)*(X-repmat(obj.xopt,size(X,1),1))*obj.R + 0.5;
            S = 100*(Z(:,1:end-1).^2-Z(:,2:end)).^2 + (Z(:,1:end-1)-1).^2;
            PopObj = 10*mean(S/4000-cos(S),2) + 10;
        end
    end
end
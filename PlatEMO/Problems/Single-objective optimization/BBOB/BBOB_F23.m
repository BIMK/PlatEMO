classdef BBOB_F23 < PROBLEM
% <2009> <single> <real> <expensive/none>
% Katsuura function
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
        Q;      % Orthogonal (rotation) matrix
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
            obj.Q = orth(rand(RandStream('mlfg6331_64','Seed',1),obj.D));
            obj.R = orth(rand(RandStream('mlfg6331_64','Seed',2),obj.D));
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            Z = (X-repmat(obj.xopt,size(X,1),1))*(obj.Q*diag(100.^(0.5*(0:obj.D-1)/(obj.D-1)))*obj.R);
            Z = repmat(Z,1,1,32);
            J = repmat(reshape(2.^(1:32),1,1,32),size(Z,1),size(Z,2),1);
            PopObj = 10/obj.D^2*prod((1+repmat(1:obj.D,size(Z,1),1).*sum(abs(J.*Z-round(J.*Z))./J,3)).^(10/obj.D^1.2),2) - 10/obj.D^2 + fpen(X);
        end
    end
end

function F = fpen(X)
    F = sum(max(0,abs(X)-5).^2,2);
end
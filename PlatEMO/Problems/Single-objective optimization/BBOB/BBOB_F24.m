classdef BBOB_F24 < PROBLEM
% <2009> <single> <real> <expensive/none>
% Lunacek bi-Rastrigin function

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
            obj.M = 1;
            if isempty(obj.D); obj.D = 30; end
            obj.xopt     = 1.25*sign(rand(RandStream('mlfg6331_64','Seed',1),1,obj.D)-0.5);
            obj.lower    = zeros(1,obj.D) - 5;
            obj.upper    = zeros(1,obj.D) + 5;
            obj.encoding = ones(1,obj.D);
            obj.Q = orth(rand(RandStream('mlfg6331_64','Seed',1),obj.D));
            obj.R = orth(rand(RandStream('mlfg6331_64','Seed',2),obj.D));
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            X = 2*sign(repmat(obj.xopt,size(X,1),1)).*X;
            Z = (X-2.5)*(obj.Q*diag(100.^(0.5*(0:obj.D-1)/(obj.D-1)))*obj.R);
            s = 1-1/(2*sqrt(obj.D+20)-8.2);
            PopObj = min(sum((X-2.5).^2,2),obj.D+s*sum((X+sqrt(5.25/s)).^2,2)) + 10*(obj.D-sum(cos(2*pi*Z),2)) + 1e4*fpen(X);
        end
    end
end

function F = fpen(X)
    F = sum(max(0,abs(X)-5).^2,2);
end
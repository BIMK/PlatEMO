classdef BBOB_F20 < PROBLEM
% <2009> <single> <real> <expensive/none>
% Schwefel function

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
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 1;
            if isempty(obj.D); obj.D = 30; end
            obj.xopt     = 4.2096874633/2*sign(rand(RandStream('mlfg6331_64','Seed',1),1,obj.D)-0.5);
            obj.lower    = zeros(1,obj.D) - 5;
            obj.upper    = zeros(1,obj.D) + 5;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            Xopt   = repmat(obj.xopt,size(X,1),1);
            X      = 2*sign(Xopt).*X;
            Z(:,1) = X(:,1);
            Z(:,2:obj.D) = X(:,2:end) + 0.25*(X(:,1:end-1)-2*abs(Xopt(:,1:end-1)));
            Z = 100*((Z-2*abs(Xopt))*diag(10.^(0.5*(0:obj.D-1)/(obj.D-1)))+2*abs(Xopt));
            PopObj = -0.01*mean(Z.*sin(sqrt(abs(Z))),2) + 4.189828872724339 + 100*fpen(Z/100);
        end
    end
end

function F = fpen(X)
    F = sum(max(0,abs(X)-5).^2,2);
end
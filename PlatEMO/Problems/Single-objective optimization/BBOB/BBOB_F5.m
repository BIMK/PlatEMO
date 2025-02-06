classdef BBOB_F5 < PROBLEM
% <2009> <single> <real> <expensive/none>
% Linear slope

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
            obj.xopt     = 5*sign(rand(RandStream('mlfg6331_64','Seed',1),1,obj.D)-0.5);
            obj.lower    = zeros(1,obj.D) - 5;
            obj.upper    = zeros(1,obj.D) + 5;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            Z    = X;
            Xopt = repmat(obj.xopt,size(X,1),1);
            Z(Xopt.*X>=25) = Xopt(Xopt.*X>=25);
            S = sign(Xopt).*10.^(repmat(0:obj.D-1,size(X,1),1)/(obj.D-1));
            PopObj = sum(5*abs(S)-S.*Z,2);
        end
    end
end
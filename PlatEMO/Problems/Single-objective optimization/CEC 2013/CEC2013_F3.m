classdef CEC2013_F3 < PROBLEM
% <single> <real> <large>
% Shifted Ackley's function

%------------------------------- Reference --------------------------------
% X. Li, K. Tang, M. N. Omidvar, Z. Yang, and K. Qin, Benchmark functions
% for the CEC'2013 special session and competition on large-scale global
% optimization, RMIT University, Australia, 2013.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        Xopt;	% Optimal decision vector
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'CEC2013.mat'),'Data');
            obj.Xopt = Data{3};
            obj.M    = 1;
            obj.D    = 1000;
            obj.lower    = zeros(1,obj.D) - 32;
            obj.upper    = zeros(1,obj.D) + 32;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = Ackley(PopDec-repmat(obj.Xopt,size(PopDec,1),1));
        end
    end
end

function F = Ackley(X)
    F = -20.*exp(-0.2.*sqrt(mean(X.^2,2))) - exp(mean(cos(2*pi*X),2)) + 20 + exp(1);
end
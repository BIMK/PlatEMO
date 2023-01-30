classdef SOP_F16 < PROBLEM
% <single> <real> <expensive/none>
% Six-hump camel-back function

%------------------------------- Reference --------------------------------
% X. Yao, Y. Liu, and G. Lin, Evolutionary programming made faster, IEEE
% Transactions on Evolutionary Computation, 1999, 3(2): 82-102.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 1;
            obj.D = 2;
            obj.lower    = zeros(1,obj.D) - 5;
            obj.upper    = zeros(1,obj.D) + 5;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = 4*PopDec(:,1).^2-2.1*PopDec(:,1).^4+PopDec(:,1).^6/3+PopDec(:,1).*PopDec(:,2)-4*PopDec(:,2).^2+4*PopDec(:,2).^4;
        end
        %% Generate the minimum objective value
        function R = GetOptimum(obj,N)
            R = -1.0316;
        end
    end
end
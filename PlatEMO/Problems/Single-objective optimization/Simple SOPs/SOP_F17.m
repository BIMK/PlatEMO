classdef SOP_F17 < PROBLEM
% <single> <real> <expensive/none>
% Branin function

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
            obj.lower    = [-5, 0];
            obj.upper    = [10,15];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = (PopDec(:,2)-5.1/4/pi.^2*PopDec(:,1).^2+5/pi*PopDec(:,1)-6).^2 + 10*(1-1/8/pi)*cos(PopDec(:,1)) + 10;
        end
        %% Generate the minimum objective value
        function R = GetOptimum(obj,N)
            R = 0.3979;
        end
    end
end
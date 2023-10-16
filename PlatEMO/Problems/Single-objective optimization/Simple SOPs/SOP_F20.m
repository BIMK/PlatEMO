classdef SOP_F20 < PROBLEM
% <single> <real> <expensive/none>
% Hartman's family

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
            obj.D = 6;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            a = [10 3 17 3.5 1.7 8;0.05 10 17 0.1 8 14;3 3.5 1.7 10 17 8;17 8 0.05 10 0.1 14];
            c = [1;1.2;3;3.2];
            p = [0.1312 0.1696 0.5569 0.0124 0.8283 0.5886;0.2329 0.4135 0.8307 0.3736 0.1004 0.9991;0.2348 0.1415 0.3522 0.2883 0.3047 0.6650;0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                PopObj(i) = -sum(c.*exp(-sum(a.*(repmat(PopDec(i,:),4,1)-p).^2,2)));
            end
        end
        %% Generate the minimum objective value
        function R = GetOptimum(obj,N)
            R = -3.322;
        end
    end
end
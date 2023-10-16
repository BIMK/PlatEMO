classdef SOP_F19 < PROBLEM
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
            obj.D = 3;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            a = [3 10 30;0.1 10 35;3 10 30;0.1 10 35];
            c = [1;1.2;3;3.2];
            p = [0.3689 0.1170 0.2673;0.4699 0.4387 0.7470;0.1091 0.8732 0.5547;0.03815 0.5743 0.8828];
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                PopObj(i) = -sum(c.*exp(-sum(a.*(repmat(PopDec(i,:),4,1)-p).^2,2)));
            end
        end
        %% Generate the minimum objective value
        function R = GetOptimum(obj,N)
            R = -3.863;
        end
    end
end
classdef SOP_F21 < PROBLEM
% <single> <real> <expensive/none>
% Shekel's family

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
            obj.D = 4;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = zeros(1,obj.D) + 10;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            a = [4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7];
            c = [0.1;0.2;0.2;0.4;0.4];
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                PopObj(i) = -sum(1./(sum((repmat(PopDec(i,:),5,1)-a).^2,2)+c));
            end
        end
        %% Generate the minimum objective value
        function R = GetOptimum(obj,N)
            R = -10.16;
        end
    end
end
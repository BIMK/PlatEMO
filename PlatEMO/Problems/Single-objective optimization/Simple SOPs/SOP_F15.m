classdef SOP_F15 < PROBLEM
% <single> <real> <expensive/none>
% Kowalik's function

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
            obj.lower    = zeros(1,obj.D) - 5;
            obj.upper    = zeros(1,obj.D) + 5;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            a = [0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
            b = 1./[0.25 0.5 1 2 4 6 8 10 12 14 16];
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                PopObj(i) = sum((a-PopDec(i,1).*(b.^2+b.*PopDec(i,2))./(b.^2+b.*PopDec(i,3)+PopDec(i,4))).^2);
            end
        end
        %% Generate the minimum objective value
        function R = GetOptimum(obj,N)
            R = 0.0003;
        end
    end
end
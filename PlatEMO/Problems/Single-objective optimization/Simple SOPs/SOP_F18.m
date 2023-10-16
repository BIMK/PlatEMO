classdef SOP_F18 < PROBLEM
% <single> <real> <expensive/none>
% Goldstein-price function

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
            obj.lower    = [-2,-2];
            obj.upper    = [ 2, 2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            PopObj = (1+(X(:,1)+X(:,2)+1).^2.*(19-14*X(:,1)+3*X(:,1).^2-14*X(:,2)+6*X(:,1).*X(:,2)+3*X(:,2).^2)).*...
                     (30+(2*X(:,1)-3*X(:,2)).^2.*(18-32*X(:,1)+12*X(:,1).^2+48*X(:,2)-36*X(:,1).*X(:,2)+27*X(:,2).^2));
        end
        %% Generate the minimum objective value
        function R = GetOptimum(obj,N)
            R = 3;
        end
    end
end
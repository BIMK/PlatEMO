classdef SOP_F14 < PROBLEM
% <single> <real> <expensive/none>
% Shekel's foxholes function

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
            obj.lower    = zeros(1,obj.D) - 65.536;
            obj.upper    = zeros(1,obj.D) + 65.536;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            a = repmat(-32:16:32,5,1);
            a = [reshape(a',1,[]);reshape(a,1,[])];
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                PopObj(i) = 1./(1/500+sum(1./((1:25)+sum((repmat(PopDec(i,:)',1,25)-a).^6,1))));
            end
        end
        %% Generate the minimum objective value
        function R = GetOptimum(obj,N)
            R = 1;
        end
    end
end
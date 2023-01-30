classdef SOP_F13 < PROBLEM
% <single> <real> <expensive/none>
% Generalized penalized function

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
            if isempty(obj.D); obj.D = 30; end
            obj.lower    = zeros(1,obj.D) - 50;
            obj.upper    = zeros(1,obj.D) + 50;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = 0.1*(sin(3*pi*PopDec(:,1)).^2+sum((PopDec(:,1:end-1)-1).^2.*(1+sin(3*pi*PopDec(:,2:end)).^2),2)+(PopDec(:,end)-1).^2.*(1+sin(2*pi*PopDec(:,end)).^2)) + sum(u(PopDec,5,100,4),2);
        end
    end
end

function X = u(X,a,k,m)
    temp     = abs(X) > a;
    X(temp)  = k.*(abs(X(temp))-a).^m;
    X(~temp) = 0;
end
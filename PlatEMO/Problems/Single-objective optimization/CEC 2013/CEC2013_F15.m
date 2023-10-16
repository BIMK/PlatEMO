classdef CEC2013_F15 < PROBLEM
% <single> <real> <large>
% Shifted Schwefel's function

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
            obj.Xopt = Data{15};
            obj.M    = 1;
            obj.D    = 1000;
            obj.lower    = zeros(1,obj.D) - 100;
            obj.upper    = zeros(1,obj.D) + 100;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = Schwefel(Tasy(Tosz(PopDec-repmat(obj.Xopt,size(PopDec,1),1)),0.2));
        end
    end
end

function F = Schwefel(X)
    F = sum(cumsum(X,2).^2,2);
end

function Z = Tosz(X)
    X1 = zeros(size(X));
    X1(X~=0) = log(abs(X(X~=0)));
    C1 = zeros(size(X)) + 5.5;
    C1(X>0) = 10;
    C2 = zeros(size(X)) + 3.1;
    C2(X>0) = 7.9;
    Z = sign(X).*exp(X1+0.049*(sin(C1.*X1)+sin(C2.*X1)));
end

function Z = Tasy(X,beta)
    Z = X.^(1+repmat(beta*linspace(0,1,size(X,2)),size(X,1),1).*sqrt(X));
    Z(X<=0) = X(X<=0);
end
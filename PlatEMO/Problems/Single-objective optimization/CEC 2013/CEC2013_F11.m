classdef CEC2013_F11 < PROBLEM
% <single> <real> <large>
% 20-nonseparable shifted and rotated Schwefel's function

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
        R25;    % Rotation matrices
        R50;
        R100;
        p;      % Rank of decision variables
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'CEC2013.mat'),'Data');
            obj.Xopt = Data{11}.xopt;
            obj.R25  = Data{11}.R25;
            obj.R50  = Data{11}.R50;
            obj.R100 = Data{11}.R100;
            obj.p    = Data{11}.p;
            obj.M    = 1;
            obj.D    = 1000;
            obj.lower    = zeros(1,obj.D) - 100;
            obj.upper    = zeros(1,obj.D) + 100;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            S = [50 50 25 25 100 100 25 25 50 25 100 25 100 50 25 25 25 100 50 25];
            W = [1.61e-2 1.29e-1 1.20e-3 3.49e-1 3.99e0 7.45e0 2.61e0 1.86e-5 7.80e-2 4.95e6 9.08e2 1.25e3 1.28e-4 2.55e-3 1.23e-2 2.25e-1 1.60e4 4.15e0 4.21e3 8.98e-6];
            PopDec = PopDec - repmat(obj.Xopt,size(PopDec,1),1);
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : length(S)
                loc = obj.p(sum(S(1:i-1))+1:sum(S(1:i)));
                switch S(i)
                    case 25
                        PopDec(:,loc) = PopDec(:,loc)*obj.R25;
                    case 50
                        PopDec(:,loc) = PopDec(:,loc)*obj.R50;
                    case 100
                        PopDec(:,loc) = PopDec(:,loc)*obj.R100;
                end
            	PopObj = PopObj + W(i)*Schwefel(Tasy(Tosz(PopDec(:,loc)),0.2));
            end
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
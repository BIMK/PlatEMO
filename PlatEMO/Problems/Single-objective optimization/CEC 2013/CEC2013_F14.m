classdef CEC2013_F14 < PROBLEM
% <single> <real> <large>
% Shifted Schwefel's function with conflicting overlapping subcomponents

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
            obj.Xopt = Data{14}.xopt;
            obj.R25  = Data{14}.R25;
            obj.R50  = Data{14}.R50;
            obj.R100 = Data{14}.R100;
            obj.p    = Data{14}.p;
            obj.M    = 1;
            obj.D    = 905;
            obj.lower    = zeros(1,obj.D) - 100;
            obj.upper    = zeros(1,obj.D) + 100;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            S = [50 50 25 25 100 100 25 25 50 25 100 25 100 50 25 25 25 100 50 25];
            W = [4.75e-1 4.99e5 3.28e2 3.23e-1 1.36e2 9.03e0 9.24e-2 1.10e-4 9.37e-3 3.00e2 4.94e0 8.14e1 6.54e-1 1.16e1 2.86e6 8.58e-5 2.36e1 4.81e-2 1.43e0 1.22e1];
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : length(S)
                loc1 = obj.p(sum(S(1:i-1))+1-(i-1)*5:sum(S(1:i))-(i-1)*5);
                loc2 = sum(S(1:i-1))+1:sum(S(1:i));
                switch S(i)
                    case 25
                        X = (PopDec(:,loc1)-repmat(obj.Xopt(:,loc2),size(PopDec,1),1))*obj.R25;
                    case 50
                        X = (PopDec(:,loc1)-repmat(obj.Xopt(:,loc2),size(PopDec,1),1))*obj.R50;
                    case 100
                        X = (PopDec(:,loc1)-repmat(obj.Xopt(:,loc2),size(PopDec,1),1))*obj.R100;
                end
            	PopObj = PopObj + W(i)*Schwefel(Tasy(Tosz(X),0.2));
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
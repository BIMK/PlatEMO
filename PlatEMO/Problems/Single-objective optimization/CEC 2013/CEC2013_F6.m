classdef CEC2013_F6 < PROBLEM
% <single> <real> <large>
% 7-nonseparable, 1-separable shifted and rotated Ackley's function

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
            obj.Xopt = Data{6}.xopt;
            obj.R25  = Data{6}.R25;
            obj.R50  = Data{6}.R50;
            obj.R100 = Data{6}.R100;
            obj.p    = Data{6}.p;
            obj.M    = 1;
            obj.D    = 1000;
            obj.lower    = zeros(1,obj.D) - 32;
            obj.upper    = zeros(1,obj.D) + 32;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            S = [50 25 25 100 50 25 25 700];
            W = [3.53e-2 5.32e-5 8.71e-1 4.95e4 8.31e-2 3.48e-5 2.82e2 1];
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
                PopObj = PopObj + W(i)*Ackley(PopDec(:,loc));
            end
        end
    end
end

function F = Ackley(X)
    F = -20.*exp(-0.2.*sqrt(mean(X.^2,2))) - exp(mean(cos(2*pi*X),2)) + 20 + exp(1);
end
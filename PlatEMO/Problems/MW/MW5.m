classdef MW5 < PROBLEM
% <problem> <MW>
% Constrained benchmark MOP proposed by Ma and Wang

%------------------------------- Reference --------------------------------
% Z. Ma and Y. Wang, Evolutionary constrained multiobjective optimization:
% Test suite construction and performance comparisons. IEEE Transactions on
% Evolutionary Computation, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function obj = MW5()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 15;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
             [N,D] = size(X);
             M     = obj.Global.M;
             g     = 1 + sum(1 - exp(-10*((X(:,M:end).^(D-M)) - 0.5 - (repmat(M:D,N,1) - 1)/(2*D)).^2),2);
             PopObj(:,1) = g.*X(:,1);
             PopObj(:,2) = g.*sqrt(1 - (PopObj(:,1)./g).^2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            l1     = atan(PopObj(:,2)./PopObj(:,1));
            l2     = 0.5*pi - 2*abs(l1-0.25*pi);
            PopCon(:,1) = PopObj(:,1).^2 + PopObj(:,2).^2 - (1.7-0.2*sin(2*l1)).^2;
            PopCon(:,2) = (1+0.5*sin(6*l2.^3)).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
            PopCon(:,3) = (1-0.45*sin(6*l2.^3)).^2 - PopObj(:,1).^2 - PopObj(:,2).^2;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = [0 1;0.3922 0.9199;0.4862 0.8739;0.5490 0.8358;0.5970 0.8023;0.6359 0.7719;0.6969 0.7174];
            P = [P;flip(P,2)];
        end
    end
end
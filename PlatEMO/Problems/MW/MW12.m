classdef MW12 < PROBLEM
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
        function obj = MW12()
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
            PopObj(:,2) = g.*(0.85 - 0.8*(PopObj(:,1)./g) - 0.08*abs(sin(3.2*pi*(PopObj(:,1)./g))));
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            PopCon(:,1) = (1 - 0.8*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2) - PopObj(:,1)/1.5))).*(1.8 - 1.125*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2)/1.8 - PopObj(:,1)/1.6)));
            PopCon(:,2) = -(1 - 0.625*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2) - PopObj(:,1)/1.6))).*(1.4 - 0.875*PopObj(:,1) - PopObj(:,2) + 0.08*sin(2*pi*(PopObj(:,2)/1.4 - PopObj(:,1)/1.6)));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 0.85 - 0.8*P(:,1) - 0.08*abs(sin(3.2*pi*P(:,1)));
            c1      = (1-0.8*P(:,1)-P(:,2)+0.08*sin(2*pi*(P(:,2)-P(:,1)/1.5))).*(1.8-1.125*P(:,1)-P(:,2)+0.08*sin(2*pi*(P(:,2)/1.8-P(:,1)/1.6)));
            invalid = c1>0;
            while any(invalid)
                P(invalid,:) = P(invalid,:).*1.001;
                c1      = (1-0.8*P(:,1)-P(:,2)+0.08*sin(2*pi*(P(:,2)-P(:,1)/1.5))).*(1.8-1.125*P(:,1)-P(:,2)+0.08*sin(2*pi*(P(:,2)/1.8-P(:,1)/1.6)));
                invalid = c1>0;
            end
        end
    end
end
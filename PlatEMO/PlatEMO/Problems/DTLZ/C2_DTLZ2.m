classdef C2_DTLZ2 < PROBLEM
% <problem> <DTLZ variant>
% Constrained DTLZ2

%------------------------------- Reference --------------------------------
% H. Jain and K. Deb, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part II:
% Handling constraints and extending to an adaptive approach, IEEE
% Transactions on Evolutionary Computation, 2014, 18(4): 602-622.
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
        function obj = C2_DTLZ2()
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = obj.Global.M + 9;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            M      = obj.Global.M;
            g      = sum((PopDec(:,M:end)-0.5).^2,2);
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            M      = obj.Global.M;
            PopObj = obj.CalObj(PopDec);
            if M == 3
                r = 0.4;
            else
                r = 0.5;
            end
            PopCon = min(min((PopObj-1).^2+repmat(sum(PopObj.^2,2),1,M)-PopObj.^2-r^2,[],2),sum((PopObj-1/sqrt(M)).^2,2)-r^2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            M = obj.Global.M;
            P = UniformPoint(N,M);
            P = P./repmat(sqrt(sum(P.^2,2)),1,M);
            if M == 3
                r = 0.4;
            else
                r = 0.5;
            end
            P(min(min((P-1).^2+repmat(sum(P.^2,2),1,M)-P.^2-r^2,[],2),sum((P-1/sqrt(M)).^2,2)-r^2)>0,:) = [];
        end
    end
end
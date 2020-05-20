classdef DC1_DTLZ3 < PROBLEM
% <problem> <DTLZ variant>
% DTLZ3 with constrains in decision space

%------------------------------- Reference --------------------------------
% K. Li, R. Chen, G. Fu, and X. Yao, Two-archive evolutionary algorithm for
% constrained multi-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2018, 23(2): 303-315.
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
        function obj = DC1_DTLZ3()
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
            [N,D]  = size(PopDec);
            M      = obj.Global.M;
            g      = 10*(D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon = 0.5 - cos(3*pi*PopDec(:,1));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            M = obj.Global.M;
            P = UniformPoint(N,M);
            c = zeros(size(P,1),M-1);
            for i = 1 : size(P,1)
                for j = 2 : M
                    temp = P(i,j)./P(i,1).*prod(c(i,M-j+2:M-1));
                    c(i,M-j+1) = sqrt(1./(1+temp^2));
                end
            end
            P(cos(3*pi*acos(c(:,1)))<0.5,:) = [];
            P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        end
    end
end
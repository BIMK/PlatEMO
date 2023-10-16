classdef DC1_DTLZ1 < PROBLEM
% <multi/many> <real> <large/none> <constrained> <expensive/none>
% DTLZ1 with constrains in decision space

%------------------------------- Reference --------------------------------
% K. Li, R. Chen, G. Fu, and X. Yao, Two-archive evolutionary algorithm for
% constrained multi-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2018, 23(2): 303-315.
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
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M + 4; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g      = 100*(obj.D-obj.M+1+sum((PopDec(:,obj.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M:end)-0.5)),2));
            PopObj = 0.5*repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(PopDec,1),1),PopDec(:,1:obj.M-1)],2)).*[ones(size(PopDec,1),1),1-PopDec(:,obj.M-1:-1:1)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon = 0.5 - cos(3*pi*PopDec(:,1));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            x = zeros(size(R,1),obj.M-1);
            for i = 1 : size(R,1)
                for j = 2 : obj.M
                    temp = R(i,j)./R(i,1).*prod(x(i,obj.M-j+2:obj.M-1));
                    x(i,obj.M-j+1) = 1./(1+temp);
                end
            end
            R(cos(3*pi*x(:,1))<0.5,:) = [];
            R = R/2;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R(:,1) = linspace(0,1,100)';
                R(:,2) = 1 - R(:,1);
                R(cos(3*pi*R(:,1))<0.5,:) = nan;
                R = R/2;
            elseif obj.M == 3
                a = linspace(0,1,10)';
                x = a*a';
                y = a*(1-a');
                z = (1-a)*ones(size(a'));
                z(cos(3*pi*a*ones(size(a')))<0.5) = nan;
                R = {x/2,y/2,z/2};
            else
                R = [];
            end
        end
    end
end
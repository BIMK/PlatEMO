classdef DC1_DTLZ3 < PROBLEM
% <multi/many> <real> <large/none> <constrained> <expensive/none>
% DTLZ3 with constrains in decision space

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
            if isempty(obj.D); obj.D = obj.M + 9; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g      = 10*(obj.D-obj.M+1+sum((PopDec(:,obj.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M:end)-0.5)),2));
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(PopDec,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(PopDec,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon = 0.5 - cos(3*pi*PopDec(:,1));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            c = zeros(size(R,1),obj.M-1);
            for i = 1 : size(R,1)
                for j = 2 : obj.M
                    temp = R(i,j)./R(i,1).*prod(c(i,obj.M-j+2:obj.M-1));
                    c(i,obj.M-j+1) = sqrt(1./(1+temp^2));
                end
            end
            R(cos(6*acos(c(:,1)))<0.5,:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                x      = linspace(0,pi/2,100)';
                R(:,1) = cos(x);
                R(:,2) = sin(x);
                R(cos(6*x)<0.5,:) = nan;
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                x = cos(a)*cos(a');
                y = cos(a)*sin(a');
                z = sin(a)*ones(size(a'));
                z(cos(6*a*ones(size(a')))<0.5) = nan;
                R = {x,y,z};
            else
                R = [];
            end
        end
    end
end
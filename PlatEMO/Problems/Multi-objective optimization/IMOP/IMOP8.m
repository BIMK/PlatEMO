classdef IMOP8 < PROBLEM
% <multi> <real> <expensive/none>
% Benchmark MOP with irregular Pareto front
% a1 --- 0.05 --- Parameter a1
% a2 ---   10 --- Parameter a2
% K  ---    5 --- Parameter K

%------------------------------- Reference --------------------------------
% Y. Tian, R. Cheng, X. Zhang, M. Li, and Y. Jin, Diversity assessment of
% multi-objective evolutionary algorithms: Performance metric and benchmark
% problems, IEEE Computational Intelligence Magazine, 2019, 14(3): 61-74.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        a1 = 0.05;  % Parameter a1
        a2 = 10;    % Parameter a2
        K  = 5;     % Parameter K
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.a1,obj.a2,obj.K] = obj.ParameterSet(0.05,10,5);
            obj.M = 3;
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            y1 = mean(PopDec(:,1:2:obj.K),2).^obj.a1;
            y2 = mean(PopDec(:,2:2:obj.K),2).^obj.a2;
            g  = sum((PopDec(:,obj.K+1:end)-0.5).^2,2);
            PopObj(:,1) = y1;
            PopObj(:,2) = y2;
            PopObj(:,3) = (1+g).*(3-sum(PopObj(:,1:2)./(1+repmat(g,1,2)).*(1+sin(19*pi.*PopObj(:,1:2))),2));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            [x,y] = meshgrid(linspace(0,1,ceil(sqrt(N))));
            R = [x(:),y(:)];
            R = [R,3-sum(R.*(1+sin(19*pi.*R)),2)];
            R = R(NDSort(R,1)==1,:);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            [x,y]  = meshgrid(linspace(0,1,150));
            z      = 3 - x.*(1+sin(19*pi.*x)) - y.*(1+sin(19*pi.*y));
            nd     = reshape(NDSort([x(:),y(:),z(:)],1)==1,size(z));
            z(~nd) = nan;
            R      = {x,y,z};
        end
    end
end
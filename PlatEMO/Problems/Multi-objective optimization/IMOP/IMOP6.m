classdef IMOP6 < PROBLEM
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
            r  = max(0,min(sin(3*pi*y1).^2,sin(3*pi*y2).^2)-0.05);
            PopObj(:,1) = (1+g).*y1 + ceil(r);
            PopObj(:,2) = (1+g).*y2 + ceil(r);
            PopObj(:,3) = (0.5+g).*(2-y1-y2) + ceil(r);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            [x,y] = meshgrid(linspace(0,1,ceil(sqrt(N))));
            R = [x(:),y(:)];
            r = max(0,min(sin(3*pi*R(:,1)).^2,sin(3*pi*R(:,2)).^2)-0.05);
            R(:,3) = 1 - sum(R,2)/2;
            R = R + repmat(ceil(r),1,3);
            R = R(NDSort(R,1)==1,:);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            [x,y] = meshgrid(linspace(0,1,50));
            z     = 1 - x/2 - y/2;
            R     = [x(:),y(:),z(:)];
            r     = max(0,min(sin(3*pi*R(:,1)).^2,sin(3*pi*R(:,2)).^2)-0.05);
            R     = R + repmat(ceil(r),1,3);
            fes   = NDSort(R,1) == 1;
            z(reshape(~fes,size(z))) = nan;
            R     = {x,y,z};
        end
    end
end
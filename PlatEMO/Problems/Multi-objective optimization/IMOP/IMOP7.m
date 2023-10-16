classdef IMOP7 < PROBLEM
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
            PopObj(:,1) = (1+g).*cos(y1*pi/2).*cos(y2*pi/2);
            PopObj(:,2) = (1+g).*cos(y1*pi/2).*sin(y2*pi/2);
            PopObj(:,3) = (1+g).*sin(y1*pi/2);
            r = min(min(abs(PopObj(:,1)-PopObj(:,2)),abs(PopObj(:,2)-PopObj(:,3))),abs(PopObj(:,3)-PopObj(:,1)));
            PopObj = PopObj + repmat(10*max(0,r-0.1),1,3);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,3);
            R = R./repmat(sqrt(sum(R.^2,2)),1,3);
            r = min(min(abs(R(:,1)-R(:,2)),abs(R(:,2)-R(:,3))),abs(R(:,3)-R(:,1)));
            R(r>0.1,:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            a = linspace(0,pi/2,50)';
            x = sin(a)*cos(a');
            y = sin(a)*sin(a');
            z = cos(a)*ones(size(a'));
            R = [x(:),y(:),z(:)];
            fes = min(min(abs(R(:,1)-R(:,2)),abs(R(:,2)-R(:,3))),abs(R(:,3)-R(:,1))) <= 0.1;
            z(reshape(~fes,size(z))) = nan;
            R = {x,y,z};
        end
    end
end
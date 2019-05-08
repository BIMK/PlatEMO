classdef IMOP5 < PROBLEM
% <problem> <IMOP>
% Benchmark MOP with irregular Pareto front
% a1 --- 0.05 --- Parameter a1
% a2 ---   10 --- Parameter a2
% K  ---    5 --- Parameter K

%------------------------------- Reference --------------------------------
% Y. Tian, R. Cheng, X. Zhang, M. Li, and Y. Jin, Diversity assessment of
% multi-objective evolutionary algorithms: Performance metric and benchmark
% problems, IEEE Computational Intelligence Magazine, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
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
        %% Initialization
        function obj = IMOP5()
            [obj.a1,obj.a2,obj.K] = obj.Global.ParameterSet(0.05,10,5);
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            y1 = mean(PopDec(:,1:2:obj.K),2).^obj.a1;
            y2 = mean(PopDec(:,2:2:obj.K),2).^obj.a2;
            g  = sum((PopDec(:,obj.K+1:end)-0.5).^2,2);
            PopObj(:,1) = 0.4*cos(pi*ceil(y1*8)/4) + 0.1*y2.*cos(16*pi*y1);
            PopObj(:,2) = 0.4*sin(pi*ceil(y1*8)/4) + 0.1*y2.*sin(16*pi*y1);
            PopObj(:,3) = 0.5 - sum(PopObj(:,1:2),2);
            PopObj = PopObj + repmat(g,1,3);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            x = ReplicatePoint(ceil(N/8*1.3),2) - 0.5;
            x = 0.2*x(sum(x.^2,2)<=0.25,:);
            r = [0.4*cos((1:8)'*pi/4),0.4*sin((1:8)'*pi/4)];
            P = [];
            for i = 1 : size(r,1)
                P = [P;x+repmat(r(i,:),size(x,1),1)];
            end
            P = [P,0.5-sum(P,2)];
        end
    end
end

function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end
classdef ZDT6 < PROBLEM
% <problem> <ZDT>
% Benchmark MOP proposed by Zitzler, Deb, and Thiele

%------------------------------- Reference --------------------------------
% E. Zitzler, K. Deb, and L. Thiele, Comparison of multiobjective
% evolutionary algorithms: Empirical results, Evolutionary computation,
% 2000, 8(2): 173-195.
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
        function obj = ZDT6()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj(:,1) = 1 - exp(-4*PopDec(:,1)).*sin(6*pi*PopDec(:,1)).^6;
            g = 1 + 9*sum(PopDec(:,2:end),2).^0.25;
            h = 1 - (PopObj(:,1)./g).^2;
            PopObj(:,2) = g.*h;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            minf1  = 0.280775;
            P(:,1) = (minf1:(1-minf1)/(N-1):1)';
            P(:,2) = 1 - P(:,1).^2;
        end
    end
end
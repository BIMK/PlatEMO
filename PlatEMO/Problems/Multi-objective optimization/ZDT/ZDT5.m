classdef ZDT5 < PROBLEM
% <multi> <binary> <large/none> <expensive/none>
% Benchmark MOP proposed by Zitzler, Deb, and Thiele

%------------------------------- Reference --------------------------------
% E. Zitzler, K. Deb, and L. Thiele, Comparison of multiobjective
% evolutionary algorithms: Empirical results, Evolutionary computation,
% 2000, 8(2): 173-195.
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
            obj.M = 2;
            if isempty(obj.D); obj.D = 80; end
            obj.D        = ceil(max(obj.D-30,1)/5)*5 + 30;
            obj.encoding = 4 + zeros(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            u      = zeros(size(PopDec,1),1+(size(PopDec,2)-30)/5);
            u(:,1) = sum(PopDec(:,1:30),2);
            for i = 2 : size(u,2)
                u(:,i) = sum(PopDec(:,(i-2)*5+31:(i-2)*5+35),2);
            end
            v           = zeros(size(u));
            v(u<5)      = 2 + u(u<5);
            v(u==5)     = 1;
            PopObj(:,1) = 1 + u(:,1);
            g           = sum(v(:,2:end),2);
            h           = 1./PopObj(:,1);
            PopObj(:,2) = g.*h;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = 1 : 31;
            R(:,2) = (obj.D-30)./5./R(:,1);
        end
    end
end
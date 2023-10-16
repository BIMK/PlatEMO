classdef FCP3 < PROBLEM
% <multi> <real> <constrained>
% Benchmark constrained MOP proposed by Jiawei Yuan

%------------------------------- Reference --------------------------------
% J. Yuan, H. Liu, Y. Ong, and Z. He, Indicator-based evolutionary
% algorithm for solving constrained multi-objective optimization problems,
% IEEE Transactions on Evolutionary Computation, 2022, 26(2): 379-391.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiawei Yuan

    methods
        %% Initialization
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D)
                obj.D = 30;
            end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g = 1 + 9*mean(PopDec(:,2:end),2);
            t = mod(floor(100*g),2);
            g = g + t.*(g-9).^2;
            PopObj(:,1) = cos(0.5*pi*PopDec(:,1)).*g;
            PopObj(:,2) = sin(0.5*pi*PopDec(:,1)).*g;
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            g = 1 + 9*mean(PopDec(:,2:end),2);
            t = mod(floor(100*g),2);
            g = g + t.*(g-9).^2;
            Dis    = abs(9-g);
            %%%%% Type-II constraints
            y1     = Dis.^2-0.25;
            y2     = 1./(Dis+1e-6).*(1.2+sin(Dis*pi));
            PopCon = min([y1,y2],[],2);
        end
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            t = 0.5*pi*(0:1/N:1)';
            P=8.5*[cos(t),sin(t)];
        end
    end
end
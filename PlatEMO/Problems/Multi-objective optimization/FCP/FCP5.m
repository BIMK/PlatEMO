classdef FCP5 < PROBLEM
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
            PopObj(:,1) = PopDec(:,1).*g;
            PopObj(:,2) = (1-PopDec(:,1)).*g;
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            g = 1 + 9*mean(PopDec(:,2:end),2); 
            %%%%% Type-III constraints
            c1     = log(sqrt((10*PopDec(:,1)-9).^2+(g-3).^2)+0.5);
            c2     = log(sqrt((10*PopDec(:,1)-6).^2+(g-6).^2)+0.05);
            c3     = (10*PopDec(:,1)-sqrt(2)).^2+(g-10).^2-2;
            c4     = 1.2+sin(pi*sqrt(c3+2));
            PopCon = min([c1,c2,c3,c4],[],2);
        end
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            t=0.5*pi*(0:1/N:1);
            c1x1 = 0.1*[9+0.5*cos(t),9-0.5*cos(t)];
            c1g  = repmat(3-0.5*sin(t),1,2);
            c2x1 = 0.1*[6+0.95*cos(t),6-0.95*cos(t)];
            c2g  = repmat(6-0.95*sin(t),1,2);
            c3x1 = 0.1*[sqrt(2)+sqrt(2)*cos(t),sqrt(2)-sqrt(2)*cos(t)];
            c3g  = repmat(10-sqrt(2)*sin(t),1,2);
            
            x1 = [c1x1,c2x1,c3x1];
            g = [c1g,c2g,c3g];
            P1 = x1.*g;
            P2 = (1-x1).*g;
            P=[P1',P2'];
            P = P(find(NDSort(P,1)==1),:);
        end
    end
end
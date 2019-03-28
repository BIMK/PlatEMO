classdef VNT2 < PROBLEM
% <problem> <VNT>
% Benchmark MOP proposed by Viennet

%------------------------------- Reference --------------------------------
% R. Viennet, C. Fonteix, and I. Marc, Multicriteria optimization using a
% genetic algorithm for determining a Pareto set, International Journal of
% Systems Science, 1996, 27(2): 255-260.
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
        function obj = VNT2()
            obj.Global.M        = 3;
            obj.Global.D        = 2;
            obj.Global.lower    = [-4,-4];
            obj.Global.upper    = [4,4];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj(:,1) = (PopDec(:,1)-2).^2/2 + (PopDec(:,2)+1).^2/13 + 3;
            PopObj(:,2) = (PopDec(:,1)+PopDec(:,2)-3).^2/36 + (-PopDec(:,1)+PopDec(:,2)+2).^2/8 - 17;
            PopObj(:,3) = (PopDec(:,1)+2*PopDec(:,2)-1).^2/175 + (2*PopDec(:,2)-PopDec(:,1)).^2/17 - 13;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            X = ReplicatePoint(N,2)*8-4;
            P = obj.CalObj(X);
            P = P(NDSort(P,1)==1,:);
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
classdef VNT4 < PROBLEM
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
        function obj = VNT4()
            obj.Global.M        = 3;
            obj.Global.D        = 2;
            obj.Global.lower    = [-4,-4];
            obj.Global.upper    = [4,4];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj(:,1) = (PopDec(:,1)-2).^2/2 + (PopDec(:,2)+1).^2/13 + 3;
            PopObj(:,2) = (PopDec(:,1)+PopDec(:,2)-3).^2/175 + (2*PopDec(:,2)-PopDec(:,1)).^2/17 - 13;
            PopObj(:,3) = (3*PopDec(:,1)-2*PopDec(:,2)+4).^2/8 + (PopDec(:,1)-PopDec(:,2)+1).^2/27 + 15;
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon(:,1) = PopDec(:,2) + 4*PopDec(:,1) - 4;
            PopCon(:,2) = -1 - PopDec(:,1);
            PopCon(:,3) = PopDec(:,1) - 2 - PopDec(:,2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            N = ceil(sqrt(N));
            x = linspace(-1,1.2,N);
            X = [];
            for i = 1 : N
                X = [X;repmat(x(i),N,1),linspace(x(i)-2,-4*x(i)+4,N)'];
            end
            P = obj.CalObj(X);
            P = P(NDSort(P,1)==1,:);
        end
    end
end
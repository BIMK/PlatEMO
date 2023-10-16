classdef VNT4 < PROBLEM
% <multi> <real> <constrained>
% Benchmark MOP proposed by Viennet

%------------------------------- Reference --------------------------------
% R. Viennet, C. Fonteix, and I. Marc, Multicriteria optimization using a
% genetic algorithm for determining a Pareto set, International Journal of
% Systems Science, 1996, 27(2): 255-260.
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
            obj.M        = 3;
            obj.D        = 2;
            obj.lower    = [-4,-4];
            obj.upper    = [4,4];
            obj.encoding = ones(1,obj.D);
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
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            N = ceil(sqrt(N));
            x = linspace(-1,1.2,N);
            X = [];
            for i = 1 : N
                X = [X;repmat(x(i),N,1),linspace(x(i)-2,-4*x(i)+4,N)'];
            end
            R = obj.CalObj(X);
            R = R(NDSort(R,1)==1,:);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            x = linspace(-1,1.2,40);
            X = [];
            for i = 1 : 40
                X = [X;repmat(x(i),40,1),linspace(x(i)-2,-4*x(i)+4,40)'];
            end
            R = obj.CalObj(X);
            R(NDSort(R,1)>1,:) = nan;
            R = {reshape(R(:,1),40,40),reshape(R(:,2),40,40),reshape(R(:,3),40,40)};
        end
    end
end
classdef VNT1 < PROBLEM
% <multi> <real>
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
            obj.lower    = [-2,-2];
            obj.upper    = [2,2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj(:,1) = PopDec(:,1).^2 + (PopDec(:,2)-1).^2;
            PopObj(:,2) = PopDec(:,1).^2 + (PopDec(:,2)+1).^2 + 1;
            PopObj(:,3) = (PopDec(:,1)-1).^2 + PopDec(:,2).^2 + 2;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            X = UniformPoint(N,2,'grid')*4-2;
            R = obj.CalObj(X);
            R = R(NDSort(R,1)==1,:);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            X = UniformPoint(1600,2,'grid')*4-2;
            R = obj.CalObj(X);
            R(NDSort(R,1)>1,:) = nan;
            R = {reshape(R(:,1),40,40),reshape(R(:,2),40,40),reshape(R(:,3),40,40)};
        end
    end
end
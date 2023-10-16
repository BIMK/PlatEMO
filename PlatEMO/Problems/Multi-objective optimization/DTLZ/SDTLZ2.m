classdef SDTLZ2 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% Scaled DTLZ2
% a --- 2 --- Scaling factor

%------------------------------- Reference --------------------------------
% K. Deb and H. Jain, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part I:
% Solving problems with box constraints, IEEE Transactions on Evolutionary
% Computation, 2014, 18(4): 577-601.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        a = 2;  % Scaling factor
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.a = obj.ParameterSet(2);
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+9; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g      = sum((PopDec(:,obj.M:end)-0.5).^2,2);
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
            PopObj = PopObj.*repmat(obj.a.^(0:obj.M-1),size(PopObj,1),1);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            R = R.*repmat(obj.a.^(0:obj.M-1),size(R,1),1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a').*obj.a,cos(a)*ones(size(a')).*obj.a.^2};
            else
                R = [];
            end
        end
    end
end
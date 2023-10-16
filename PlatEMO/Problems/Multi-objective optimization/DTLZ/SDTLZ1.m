classdef SDTLZ1 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% Scaled DTLZ1
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
        a = 2;	% Scaling factor
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.a = obj.ParameterSet(2);
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+4; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g      = 100*(obj.D-obj.M+1+sum((PopDec(:,obj.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M:end)-0.5)),2));
            PopObj = 0.5*repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(PopDec,1),1),PopDec(:,1:obj.M-1)],2)).*[ones(size(PopDec,1),1),1-PopDec(:,obj.M-1:-1:1)];
            PopObj = PopObj.*repmat(obj.a.^(0:obj.M-1),size(PopObj,1),1);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M)/2;
            R = R.*repmat(obj.a.^(0:obj.M-1),size(R,1),1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,1,10)';
                R = {a*a'/2,a*(1-a')/2.*obj.a,(1-a)*ones(size(a'))/2.*obj.a^2};
            else
                R = [];
            end
        end
    end
end
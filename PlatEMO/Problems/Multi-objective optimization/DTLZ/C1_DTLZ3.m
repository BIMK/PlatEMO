classdef C1_DTLZ3 < PROBLEM
% <multi/many> <real> <large/none> <constrained> <expensive/none>
% Constrained DTLZ3

%------------------------------- Reference --------------------------------
% H. Jain and K. Deb, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part II:
% Handling constraints and extending to an adaptive approach, IEEE
% Transactions on Evolutionary Computation, 2014, 18(4): 602-622.
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
            if isempty(obj.M)
                obj.M = 3;
            end
            if isempty(obj.D)
                obj.D = obj.M + 9;
            end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            PopDec = varargin{1};
            PopDec = max(min(PopDec,repmat(obj.upper,size(PopDec,1),1)),repmat(obj.lower,size(PopDec,1),1));
            g      = 10*(obj.D-obj.M+1+sum((PopDec(:,obj.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M:end)-0.5)),2));
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(PopDec,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(PopDec,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
            if obj.M == 2
                r = 6;
            elseif obj.M <= 3
                r = 9;
            elseif obj.M <= 8
                r = 12.5;
            else
                r = 15;
            end
            PopCon = -(sum(PopObj.^2,2)-16).*(sum(PopObj.^2,2)-r.^2);
            Population = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
    end
end
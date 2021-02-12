classdef C3_DTLZ4 < PROBLEM
% <multi/many> <real> <large/none> <constrained> <expensive/none>
% Constrained DTLZ4

%------------------------------- Reference --------------------------------
% H. Jain and K. Deb, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part II:
% Handling constraints and extending to an adaptive approach, IEEE
% Transactions on Evolutionary Computation, 2014, 18(4): 602-622.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
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
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec(:,1:obj.M-1) = PopDec(:,1:obj.M-1).^100;
            g      = sum((PopDec(:,obj.M:end)-0.5).^2,2);
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopObj = obj.CalObj(PopDec);
            PopCon = 1 - PopObj.^2/4 - (repmat(sum(PopObj.^2,2),1,obj.M)-PopObj.^2);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)-3/4*max(R.^2,[],2)),1,obj.M);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                x = sin(a)*cos(a');
                y = sin(a)*sin(a');
                z = cos(a)*ones(size(a'));
                R = [x(:),y(:),z(:)];
                R = R./repmat(sqrt(sum(R.^2,2)-3/4*max(R.^2,[],2)),1,3);
                R = {reshape(R(:,1),size(x)),reshape(R(:,2),size(x)),reshape(R(:,3),size(x))};
            else
                R = [];
            end
        end
    end
end
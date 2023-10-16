classdef C2_DTLZ2 < PROBLEM
% <multi/many> <real> <large/none> <constrained> <expensive/none>
% Constrained DTLZ2

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
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M + 9; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            PopDec = varargin{1};
            PopDec = max(min(PopDec,repmat(obj.upper,size(PopDec,1),1)),repmat(obj.lower,size(PopDec,1),1));
            g      = sum((PopDec(:,obj.M:end)-0.5).^2,2);
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
            if obj.M == 3
                r = 0.4;
            else
                r = 0.5;
            end
            PopCon = min(min((PopObj-1).^2+repmat(sum(PopObj.^2,2),1,obj.M)-PopObj.^2-r^2,[],2),sum((PopObj-1/sqrt(obj.M)).^2,2)-r^2);
            Population = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            if obj.M == 3
                r = 0.4;
            else
                r = 0.5;
            end
            R(min(min((R-1).^2+repmat(sum(R.^2,2),1,obj.M)-R.^2-r^2,[],2),sum((R-1/sqrt(obj.M)).^2,2)-r^2)>0,:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
                R(min(min((R-1).^2+repmat(sum(R.^2,2),1,2)-R.^2-0.5^2,[],2),sum((R-1/sqrt(2)).^2,2)-0.5^2)>0,:) = nan;
            elseif obj.M == 3
                a = linspace(0,pi/2,30)';
                x = sin(a)*cos(a');
                y = sin(a)*sin(a');
                z = cos(a)*ones(size(a'));
                R = [x(:),y(:),z(:)];
                fes = min(min((R-1).^2+repmat(sum(R.^2,2),1,3)-R.^2-0.4^2,[],2),sum((R-1/sqrt(3)).^2,2)-0.4^2) <= 0;
                z(reshape(~fes,size(z))) = nan;
                R = {x,y,z};
            else
                R = [];
            end
        end
    end
end
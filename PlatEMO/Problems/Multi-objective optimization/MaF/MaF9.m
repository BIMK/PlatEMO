classdef MaF9 < PROBLEM
% <multi/many> <real>
% ML-DMP

%------------------------------- Reference --------------------------------
% R. Cheng, M. Li, Y. Tian, X. Zhang, S. Yang, Y. Jin, and X. Yao, A
% benchmark test suite for evolutionary many-objective optimization,
% Complex & Intelligent Systems, 2017, 3(1): 67-81.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Points;     % Vertexes
        Polygons;   % Infeasible polygons
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            if isempty(obj.M); obj.M = 10; end
            obj.M        = max(obj.M,3);
            obj.D        = 2;
            obj.lower    = [-10000,-10000];
            obj.upper    = [10000,10000];
            obj.encoding = ones(1,obj.D);
            % Generate vertexes
            obj.Points = [];
            [thera,rho] = cart2pol(0,1);
            [obj.Points(:,1),obj.Points(:,2)] = pol2cart(thera-(1:obj.M)*2*pi/obj.M,rho);
            % Generate invalid polygons
            head = repmat((1:obj.M)',ceil(obj.M/2-2),1);
            tail = repmat(1:ceil(obj.M/2-2),obj.M,1);
            tail = head + tail(:);
            obj.Polygons = cell(1,length(head));
            for i = 1 : length(obj.Polygons)
                obj.Polygons{i} = obj.Points(mod((head(i):tail(i))-1,obj.M)+1,:);
                obj.Polygons{i} = [obj.Polygons{i};repmat(2*Intersection(obj.Points(mod([head(i)-1,head(i),tail(i),tail(i)+1]-1,obj.M)+1,:)),size(obj.Polygons{i},1),1)-obj.Polygons{i}];
            end
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj,PopDec)
            Invalid = getInvalid(PopDec,obj.Polygons,obj.Points);
            while any(Invalid)
                PopDec(Invalid,:) = unifrnd(repmat(obj.lower,sum(Invalid),1),repmat(obj.upper,sum(Invalid),1));
                Invalid           = getInvalid(PopDec,obj.Polygons,obj.Points);
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),size(obj.Points,1));
            for m = 1 : size(obj.Points,1)
                PopObj(:,m) = Point2Line(PopDec,obj.Points(mod(m-1:m,size(obj.Points,1))+1,:));
            end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            [X,Y] = ndgrid(linspace(-1,1,ceil(sqrt(N))));
            ND    = inpolygon(X(:),Y(:),obj.Points(:,1),obj.Points(:,2));
            R     = obj.CalObj([X(ND),Y(ND)]);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 3
                [X,Y]    = ndgrid(linspace(-1,1,40));
                R        = obj.CalObj([X(:),Y(:)]);
                ND       = inpolygon(X(:),Y(:),obj.Points(:,1),obj.Points(:,2));
                R(~ND,:) = nan;
                R = {reshape(R(:,1),size(X)),reshape(R(:,2),size(X)),reshape(R(:,3),size(X))};
            else
                R = [];
            end
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            Draw(obj.Points([1:end,1],:),'-k','LineWidth',1.5,{'\it x\rm_1','\it x\rm_2',[]});
            Draw(Population.decs);
        end
    end
end

function r = Intersection(p)
    if p(1,1) == p(2,1)
        r(1) = p(1,1);
        r(2) = p(3,2)+(r(1)-p(3,1))*(p(3,2)-p(4,2))/(p(3,1)-p(4,1));
    elseif p(3,1) == p(4,1)
        r(1) = p(3,1);
        r(2) = p(1,2)+(r(1)-p(1,1))*(p(1,2)-p(2,2))/(p(1,1)-p(2,1));
    else
        k1   = (p(1,2)-p(2,2))/(p(1,1)-p(2,1));
        k2   = (p(3,2)-p(4,2))/(p(3,1)-p(4,1));
        r(1) = (k1*p(1,1)-k2*p(3,1)+p(3,2)-p(1,2))/(k1-k2);
        r(2) = p(1,2)+(r(1)-p(1,1))*k1;
    end
end

function Invalid = getInvalid(PopDec,Polygons,Points)
    Invalid = false(size(PopDec,1),1);
    for i = 1 : length(Polygons)
        Invalid = Invalid | inpolygon(PopDec(:,1),PopDec(:,2),Polygons{i}(:,1),Polygons{i}(:,2));
    end
    Invalid = Invalid & ~inpolygon(PopDec(:,1),PopDec(:,2),Points(:,1),Points(:,2));
end

function Distance = Point2Line(PopDec,Line)
    Distance = abs((Line(1,1)-PopDec(:,1)).*(Line(2,2)-PopDec(:,2))-(Line(2,1)-PopDec(:,1)).*(Line(1,2)-PopDec(:,2)))./sqrt((Line(1,1)-Line(2,1)).^2+(Line(1,2)-Line(2,2)).^2);
end
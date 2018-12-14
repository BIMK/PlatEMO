classdef MLDMP < PROBLEM
% <problem> <Distance minimization>
% The multi-line distance minimization problem
% lower --- -100 --- Lower bound of decision variables
% upper ---  100 --- Upper bound of decision variables

%------------------------------- Reference --------------------------------
% M. Li, C. Grosan, S. Yang, X. Liu, and X. Yao, Multiline distance
% minimization: A visualized many-objective test problem suite, IEEE
% Transactions on Evolutionary Computation, 2018, 22(1): 61-78.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
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
        %% Initialization
        function obj = MLDMP()
            % Parameter setting
            [lower,upper] = obj.Global.ParameterSet(-100,100);
            if isempty(obj.Global.M)
                obj.Global.M = 10;
            end
            obj.Global.D        = 2;
            obj.Global.lower    = zeros(1,2) + lower;
            obj.Global.upper    = zeros(1,2) + upper;
            obj.Global.encoding = 'real';
            % Generate vertexes
            if mod(obj.Global.M,2) == 0
                Angle = (2.*(1:obj.Global.M)-3).*pi./obj.Global.M;
            else
                Angle = (2.*(1:obj.Global.M)-2).*pi./obj.Global.M;
            end
            obj.Points = [sin(Angle)',cos(Angle)'];
            % Generate infeasible polygons
            head = repmat((1:obj.Global.M)',ceil(obj.Global.M/2-2),1);
            tail = repmat(1:ceil(obj.Global.M/2-2),obj.Global.M,1);
            tail = head + tail(:);
            obj.Polygons = cell(1,length(head));
            for i = 1 : length(obj.Polygons)
                obj.Polygons{i} = obj.Points(mod((head(i):tail(i))-1,obj.Global.M)+1,:);
                obj.Polygons{i} = [obj.Polygons{i};repmat(2*Intersection(obj.Points(mod([head(i)-1,head(i),tail(i),tail(i)+1]-1,obj.Global.M)+1,:)),size(obj.Polygons{i},1),1)-obj.Polygons{i}];
            end
        end
        %% Repair infeasible solutions
        function PopDec = CalDec(obj,PopDec)
            Infeasible = getInfeasible(PopDec,obj.Polygons,obj.Points);
            while any(Infeasible)
                PopDec(Infeasible,:) = unifrnd(repmat(obj.Global.lower,sum(Infeasible),1),repmat(obj.Global.upper,sum(Infeasible),1));
                Infeasible           = getInfeasible(PopDec,obj.Polygons,obj.Points);
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),size(obj.Points,1));
            for m = 1 : size(obj.Points,1)
                PopObj(:,m) = Point2Line(PopDec,obj.Points(mod(m-1:m,size(obj.Points,1))+1,:));
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            [X,Y] = ndgrid(linspace(-1,1,ceil(sqrt(N))));
            ND    = inpolygon(X(:),Y(:),obj.Points(:,1),obj.Points(:,2));
            P     = obj.CalObj([X(ND),Y(ND)]);
        end
        %% Draw special figure
        function Draw(obj,PopDec)
            cla; Draw(PopDec);
            plot(obj.Points([1:end,1],1),obj.Points([1:end,1],2),'-k','LineWidth',1.5);
            xlabel('\itx\rm_1'); ylabel('\itx\rm_2');
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

function Infeasible = getInfeasible(PopDec,Polygons,Points)
    Infeasible = false(size(PopDec,1),1);
    for i = 1 : length(Polygons)
        Infeasible = Infeasible | inpolygon(PopDec(:,1),PopDec(:,2),Polygons{i}(:,1),Polygons{i}(:,2));
    end
    Infeasible = Infeasible & ~inpolygon(PopDec(:,1),PopDec(:,2),Points(:,1),Points(:,2));
end

function Distance = Point2Line(PopDec,Line)
    Distance = abs((Line(1,1)-PopDec(:,1)).*(Line(2,2)-PopDec(:,2))-(Line(2,1)-PopDec(:,1)).*(Line(1,2)-PopDec(:,2)))./sqrt((Line(1,1)-Line(2,1)).^2+(Line(1,2)-Line(2,2)).^2);
end
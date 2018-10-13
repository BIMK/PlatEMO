function varargout = MaF9(Operation,Global,input)
% <problem> <MaF>
% A benchmark test suite for evolutionary many-objective optimization
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

persistent Points Polygons;

    % This problem is multi-line distance minimization problem
    switch Operation
        case 'init'
            if isempty(Global.lower)
                Global.M          = 10;
                Global.D          = 2;
                Global.D          = 2;
                Global.lower      = [-10000,-10000];
                Global.upper      = [10000,10000];
                Global.operator   = @EAreal;
                Global.evaluation = max(1e5,1e4*Global.D);
                % Feasible polygon
                Points = [];
                [thera,rho] = cart2pol(0,1);
                [Points(:,1),Points(:,2)] = pol2cart(thera-(1:Global.M)*2*pi/Global.M,rho);
                % Infeasible polygons
                head     = repmat((1:Global.M)',ceil(Global.M/2-2),1);
                tail     = repmat(1:ceil(Global.M/2-2),Global.M,1);
                tail     = head + tail(:);
                Polygons = cell(1,length(head));
                for i = 1 : length(Polygons)
                    Polygons{i} = Points(mod((head(i):tail(i))-1,Global.M)+1,:);
                    Polygons{i} = [Polygons{i};repmat(2*Intersection(Points(mod([head(i)-1,head(i),tail(i),tail(i)+1]-1,Global.M)+1,:)),size(Polygons{i},1),1)-Polygons{i}];
                end
            end
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            PopDec    = MaF9('value',Global,PopDec);
            varargout = {PopDec};
        case 'value'
            PopDec     = input;
            Infeasible = getInfeasible(PopDec,Polygons,Points);
            while any(Infeasible)
                PopDec(Infeasible,:) = rand(sum(Infeasible),Global.D).*repmat(Global.upper-Global.lower,sum(Infeasible),1) + repmat(Global.lower,sum(Infeasible),1);
                Infeasible           = getInfeasible(PopDec,Polygons,Points);
            end
            
            PopObj = zeros(size(PopDec,1),size(Points,1));
            for m = 1 : size(Points,1)
                PopObj(:,m) = Point2Line(PopDec,Points(mod(m-1:m,size(Points,1))+1,:));
            end
            
            PopCon = [];
            
            varargout = {PopDec,PopObj,PopCon};
        case 'PF'
            [X,Y]      = ndgrid(linspace(-1,1,ceil(sqrt(input))));
            ND         = inpolygon(X(:),Y(:),Points(:,1),Points(:,2));
            [~,PopObj] = MaF9('value',Global,[X(ND),Y(ND)]);
            varargout  = {PopObj};
        case 'draw'
            cla; Draw(input);
            plot(Points([1:end,1],1),Points([1:end,1],2),'-k','LineWidth',1.5);
            xlabel('\itx\rm_1'); ylabel('\itx\rm_2');
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
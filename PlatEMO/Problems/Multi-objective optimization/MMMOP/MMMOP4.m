classdef MMMOP4 < PROBLEM
% <multi/many> <real> <multimodal>
% Multi-modal multi-objective optimization problem
% kA --- 1 --- Number of decision variables in XA
% c  --- 2 --- Parameter c
% d  --- 2 --- Parameter d

%------------------------------- Reference --------------------------------
% Y. Liu, G. G. Yen, and D. Gong, A multi-modal multi-objective
% evolutionary algorithm using two-archive and recombination strategies,
% IEEE Transactions on Evolutionary Computation, 2019, 23(4): 660-674.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        kA;     % Number of decision variables in XA
        c;      % Parameter c
        d;      % Parameter d
        POS;    % Pareto optimal set for IGDX calculation
    end
	methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = obj.M + 1; end
            [obj.kA,obj.c,obj.d] = obj.ParameterSet(1,2,2);
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            y    = PopDec(:,1:obj.M-1)*obj.d*(obj.d+1)/2;
            temp = cumsum(obj.d:-1:1);
            for i = 1 : size(y,1)
                for j = 1 : size(y,2)
                    index  = find(y(i,j)<=temp,1);
                    y(i,j) = (y(i,j)-temp(index)+obj.d-index+1)/obj.d;
                end
            end
            g      = 100*(obj.D-obj.M+1+sum(cos(2.*pi.*obj.c.*PopDec(:,obj.M:obj.M-1+obj.kA)),2)+sum((PopDec(:,obj.M+obj.kA:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M+obj.kA:end)-0.5)),2));          
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(y*pi/2)],2)).*[ones(size(g,1),1),sin(y(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj,N)
            % Generate points in Pareto optimal set
            XA      = Grid(0.5/obj.c:1/obj.c:(1-0.5/obj.c),obj.kA);
            X       = UniformPoint(N/size(XA,1),obj.M-1,'grid');
            obj.POS = [repmat(X,size(XA,1),1),XA(repmat(1:end,size(X,1),1),:),zeros(size(X,1)*size(XA,1),obj.D-size(X,2)-size(XA,2))+0.5];
            % Generate points on Pareto front
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                a = linspace(0,pi/2,100)';
                R = [sin(a),cos(a)];
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            switch metName
                case 'IGDX'
                    score = feval(metName,Population,obj.POS);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            PopDec = Population.decs;
            XA = Grid(0.5/obj.c:1/obj.c:(1-0.5/obj.c),obj.kA);
            a  = cumsum(obj.d:-1:1)/sum(obj.d:-1:1);
            X(obj.d) = a(obj.d);
            for i = obj.d-1 : -1 : 1
                X(i) = 2*a(i) - X(i+1);
            end
            X  = Grid(X,obj.M-1);
            R  = [repmat(X,size(XA,1),1),XA(repmat(1:end,size(X,1),1),:)];
            [~,Label]  = min(pdist2(PopDec(:,1:obj.M+obj.kA-1),R),[],2);
            tempStream = RandStream('mlfg6331_64','Seed',2);
            if obj.M == 2
                for i = 1 : size(R,1)
                    color = rand(tempStream,1,3);
                    Draw(Population(Label==i).objs+(i-1)*0.05,'o','MarkerSize',6,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                    Draw(obj.PF+(i-1)*0.05,'-','LineWidth',1,'Color',color);
                end
            elseif obj.M == 3
                for i = 1 : size(R,1)
                    color = rand(tempStream,1,3);
                    ax = Draw(Population(Label==i).objs+(i-1)*0.05,'o','MarkerSize',8,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
                    surf(ax,obj.PF{1}+(i-1)*0.05,obj.PF{2}+(i-1)*0.05,obj.PF{3}+(i-1)*0.05,'EdgeColor',color,'FaceColor','none');
                end
            else
                for i = 1 : size(R,1)
                    Draw(Population(Label==i).objs,'-','Color',rand(tempStream,1,3),'LineWidth',2);
                end
            end
        end
    end
end

function W = Grid(gap,M)
    if M < 1
        W = zeros(1,0);
    else
        eval(sprintf('[%s]=ndgrid(gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    end
end
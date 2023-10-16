classdef MMMOP6 < PROBLEM
% <multi/many> <real> <multimodal>
% Multi-modal multi-objective optimization problem
% kA --- 2 --- Number of decision variables in XA
% c  --- 2 --- Parameter c

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
            if isempty(obj.D); obj.D = obj.M + 2; end
            [obj.kA,obj.c] = obj.ParameterSet(2,2);
            obj.kA       = ceil(obj.kA/2)*2;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            y = 12*(PopDec(:,obj.M:obj.M+obj.kA-1)-0.5);
            z = 2*obj.c*PopDec(:,obj.M+obj.kA:end) - 2*floor(obj.c*PopDec(:,obj.M+obj.kA:end)) - 1;
            t = prod(sin(2*pi*repmat(reshape(PopDec(:,1:obj.M-1),size(PopDec,1),1,[]),1,obj.D-obj.M-obj.kA+1)+repmat(reshape((0:obj.D-obj.M-obj.kA)*pi/(obj.D-obj.M-obj.kA+1),1,[],1),size(PopDec,1),1,obj.M-1)),3);
            g = sum((y(:,1:2:obj.kA-1).^2+y(:,2:2:obj.kA)-11).^2+(y(:,1:2:obj.kA-1)+y(:,2:2:obj.kA).^2-7).^2,2) + sum((z-t).^2,2);
            PopObj = repmat(g,1,obj.M) + fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj,N)
            % Generate points in Pareto optimal set
            XA = [Grid([0.7500 0.2662 0.1851 0.7987],obj.kA/2),Grid([0.6667 0.7609 0.2264 0.3460],obj.kA/2)];
            XA = XA(:,[1:2:end,2:2:end]);
            XB = Grid(1:obj.c,obj.D-obj.M-obj.kA+1);
            X  = UniformPoint(N/size(XA,1)/size(XB,1),obj.M-1,'grid');
            R  = [repmat(X,size(XA,1),1),XA(repmat(1:end,size(X,1),1),:)];
            R  = [repmat(R,size(XB,1),1),XB(repmat(1:end,size(R,1),1),:)];
            t  = prod(sin(2*pi*repmat(reshape(R(:,1:obj.M-1),size(R,1),1,[]),1,obj.D-obj.M-obj.kA+1)+repmat(reshape((0:obj.D-obj.M-obj.kA)*pi/(obj.D-obj.M-obj.kA+1),1,[],1),size(R,1),1,obj.M-1)),3);
            R(:,obj.M+obj.kA:end) = (t+1)/2/obj.c + (R(:,obj.M+obj.kA:end)-1)/obj.c;
            obj.POS = R;
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
            XA = [Grid([0.7500 0.2662 0.1851 0.7987],obj.kA/2),Grid([0.6667 0.7609 0.2264 0.3460],obj.kA/2)];
            XA = XA(:,[1:2:end,2:2:end]);
            XB = Grid(1:obj.c,obj.D-obj.M-obj.kA+1);
            R  = [repmat(XA,size(XB,1),1),XB(repmat(1:end,size(XA,1),1),:)];
            t  = prod(sin(2*pi*repmat(reshape(PopDec(:,1:obj.M-1),size(PopDec,1),1,[]),1,obj.D-obj.M-obj.kA+1)+repmat(reshape((0:obj.D-obj.M-obj.kA)*pi/(obj.D-obj.M-obj.kA+1),1,[],1),size(PopDec,1),1,obj.M-1)),3);
            PopDec(:,obj.M+obj.kA:end) = (PopDec(:,obj.M+obj.kA:end)-(t+1)/2/obj.c)*obj.c + 1;
            [~,Label]  = min(pdist2(PopDec(:,obj.M:end),R),[],2);
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
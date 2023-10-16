classdef SMMOP5 < PROBLEM
% <multi/many> <real> <large/none> <multimodal> <sparse/none>
% Sparse multi-modal multi-objective optimization problem
% theta --- 0.1 --- Sparsity of the Pareto sets
% np    ---   4 --- Number of the Pareto sets

%------------------------------- Reference --------------------------------
% Y. Tian, R. Liu, X. Zhang, H. Ma, K. C. Tan, and Y. Jin, A
% multipopulation evolutionary algorithm for solving large-scale multimodal
% multiobjective optimization problems, IEEE Transactions on Evolutionary
% Computation, 2021, 25(3): 405-418.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        theta = 0.1;    % Sparsity of the Pareto sets
        np    = 4;    	% Number of the Pareto sets
        POS;            % Pareto optimal set for IGDX calculation
    end 
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.theta,obj.np] = obj.ParameterSet(0.1,4);
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 100; end
            obj.lower    = [zeros(1,obj.M-1)+0,zeros(1,obj.D-obj.M+1)-1];
            obj.upper    = [zeros(1,obj.M-1)+1,zeros(1,obj.D-obj.M+1)+2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,D] = size(X);
            M     = obj.M;   
            S     = ceil(obj.theta*(D-M));
            g     = zeros(N,obj.np);
            for i = 1 : obj.np
                g(:,i) = sum(g1(X(:,M+(i-1)*S:M+i*S-1),pi/3),2)+sum(g5(X(:,[M:M+(i-1)*S-1,M+i*S:end]),0),2);
            end
            PopObj = repmat(1+min(g,[],2)/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),1-cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),1-sin(X(:,M-1:-1:1)*pi/2)];
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj,N)
            % Generate points in Pareto optimal set
            A       = GetPS(obj.D-obj.M+1,obj.np,ceil(obj.theta*(obj.D-obj.M)));
            X       = UniformPoint(N/size(A,1),obj.M-1,'grid');
            obj.POS = [repmat(X,size(A,1),1),A(repmat(1:end,size(X,1),1),:)];
            % Generate points on Pareto front
            R = UniformPoint(N,obj.M);
            c = ones(size(R,1),obj.M);
            for i = 1 : size(R,1) 
                for j = 2 : obj.M
                    temp = R(i,j)/R(i,1)*prod(1-c(i,obj.M-j+2:obj.M-1));
                    c(i,obj.M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                end
            end
            x = acos(c)*2/pi;
            R = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                a = linspace(0,pi/2,100)';
                R = [1-cos(a),1-sin(a)];
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {(1-cos(a))*(1-cos(a')),(1-cos(a))*(1-sin(a')),(1-sin(a))*ones(size(a'))};
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
            A      = GetPS(obj.D-obj.M+1,obj.np,ceil(obj.theta*(obj.D-obj.M)));
            [~,Label]  = min(pdist2(PopDec(:,obj.M:end),A),[],2);
            tempStream = RandStream('mlfg6331_64','Seed',2);
            if obj.M == 2
                for i = 1 : size(A,1)
                    color = rand(tempStream,1,3);
                    Draw(Population(Label==i).objs+(i-1)*0.05,'o','MarkerSize',6,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                    Draw(obj.PF+(i-1)*0.05,'-','LineWidth',1,'Color',color);
                end
            elseif obj.M == 3
                for i = 1 : size(A,1)
                    color = rand(tempStream,1,3);
                    ax = Draw(Population(Label==i).objs+(i-1)*0.05,'o','MarkerSize',8,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
                    surf(ax,obj.PF{1}+(i-1)*0.05,obj.PF{2}+(i-1)*0.05,obj.PF{3}+(i-1)*0.05,'EdgeColor',color,'FaceColor','none');
                end
            else
                for i = 1 : size(A,1)
                    Draw(Population(Label==i).objs,'-','Color',rand(tempStream,1,3),'LineWidth',2);
                end
            end
        end
    end
end

function g = g1(x,t)
    g = (x-t).^2;
end

function g = g5(x,t)
    g = exp(log(2)*((x-t).^2)).*(sin(6*pi*(x-t)).^2)+(x-t).^2;
end

function PS = GetPS(D,np,S)
    PS = zeros(np,D);
    for i = 1 : np
        PS(i,(i-1)*S+1:i*S) = 1;
    end
end
classdef DSMOP12 < PROBLEM
% <2025> <multi/many> <real> <large/none> <sparse> <dynamic>
% Dynamic sparse multi-objective optimization problem
% theta --- 0.1 --- Sparsity of the Pareto set
% taut  ---  10 --- Number of generations for static optimization
% nt    ---  10 --- Number of distinct steps

%------------------------------- Reference --------------------------------
% P. Zhang, R. Zhang, Y. Tian, K. C. Tan, and X. Zhang. A dual model-based
% evolutionary framework for dynamic large-scale sparse multiobjective
% optimization. Swarm and Evolutionary Computation, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        theta = 0.1;    % Sparsity of the Pareto set
        taut;           % Number of generations for static optimization
        nt;             % Number of distinct steps 
        Optimums;       % Point sets on all Pareto fronts
        last = 0;       % Number of changed environments
        K;              % Sparsity
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.theta,obj.nt,obj.taut] = obj.ParameterSet(0.1,10,10);
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 100; end
            obj.lower    = [zeros(1,obj.M-1)+0,zeros(1,obj.D-obj.M+1)-1];
            obj.upper    = [zeros(1,obj.M-1)+1,zeros(1,obj.D-obj.M+1)+2];
            obj.K        = ceil(obj.theta*(obj.D-obj.M+1));
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate solutions
        function Population = Evaluation(obj,varargin)
            PopDec     = obj.CalDec(varargin{1});
            PopObj     = obj.CalObj(PopDec);
            PopCon     = obj.CalCon(PopDec);
            % Attach the current number of function evaluations to solutions
            Population = SOLUTION(PopDec,PopObj,PopCon,zeros(size(PopDec,1),1)+obj.FE);
            obj.FE     = obj.FE + length(Population);
        end        
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            T = floor(obj.FE/obj.N/obj.taut);
            t = T/obj.nt;
            if T > obj.last
                obj.last = T;
                obj.K    = ceil((rand*(obj.D-obj.M+1))* obj.theta);
            end        
            g        = g4(X(:,obj.M:end),repmat(linspace(0,1,obj.D-obj.M+1),size(X,1),1));
            [g,rank] = sort(g,2);
            temp     = false(size(rank));
            for i = 1 : size(rank,1)
                temp(i,X(i,obj.M-1+rank(i,:))==0) = true;
            end
            temp(:,1:obj.K) = false;
            g(temp) = 0;
            g       = sum(g,2);
            H       = fliplr(cumprod([ones(size(X,1),1),cos(X(:,1:obj.M-1)*pi/2)],2)).*[ones(size(X,1),1),sin(X(:,obj.M-1:-1:1)*pi/2)];
            PopObj  = repmat(1+g/(obj.D-obj.M+1),1,obj.M).*(H.^(1.25-0.25*sin(0.5.*pi.*t)));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            % Generate point sets on all Pareto fronts
            t = floor(0:obj.maxFE/obj.N/obj.taut)/obj.nt;
            G = 1.25-0.25*sin(0.5.*pi.*t);
            G = unique(round(G*1e6)/1e6);  
            R = UniformPoint(N,obj.M);
            c = ones(size(R,1),obj.M);            
            for i = 1 : size(R,1) 
                for j = 2 : obj.M
                    temp = R(i,j)/R(i,1)*prod(1-c(i,obj.M-j+2:obj.M-1));
                    c(i,obj.M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                end
            end
            x = acos(c)*2/pi;
            obj.Optimums = {};
            for i = 1 : length(G)
                g = G(i);
                h = (fliplr(cumprod([ones(size(x,1),1),cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),sin(x(:,end-1:-1:1)*pi/2)]);
                obj.Optimums(i,:) = {G(i),h.^g};
            end
            % Combine all point sets
            R = cat(1,obj.Optimums{:,2});
        end
        function score = CalMetric(obj,metName,Population)
            t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
            H      = 1.25-0.25*sin(0.5.*pi.*t);
            H      = round(H*1e6)/1e6;
            change = [0;find(H(1:end-1)~=H(2:end));length(H)];
            Scores = zeros(1,length(change)-1);
            allH   = cell2mat(obj.Optimums(:,1));
            for i = 1 : length(change)-1
                subPop    = Population(change(i)+1:change(i+1));
                Scores(i) = feval(metName,subPop,obj.Optimums{find(H(change(i)+1)==allH,1),2});
            end
            score = mean(Scores);
        end 
        function DrawObj(obj,Population)
            t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
            H      = 1.25-0.25*sin(0.5.*pi.*t);
            H      = round(H*1e6)/1e6;
            change = [0;find(H(1:end-1)~=H(2:end));length(H)];
            allH   = cell2mat(obj.Optimums(:,1));
            tempStream = RandStream('mlfg6331_64','Seed',2);
            for i = 1 : length(change)-1
                color = rand(tempStream,1,3);
                Draw(Population(change(i)+1:change(i+1)).objs,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                Draw(obj.Optimums{find(H(change(i)+1)==allH,1),2},'-','LineWidth',1,'Color',color);
            end
        end   
    end
end

function g = g4(x,t)
    g = (x-pi/3).^2 + t.*sin(6*pi*(x-pi/3)).^2;
end
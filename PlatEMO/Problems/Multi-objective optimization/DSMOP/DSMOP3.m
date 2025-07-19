classdef DSMOP3 < PROBLEM
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
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.theta,obj.nt,obj.taut] = obj.ParameterSet(0.1,10,10);
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 100; end
            obj.lower    = [zeros(1,obj.M-1)+0,zeros(1,obj.D-obj.M+1)-1];
            obj.upper    = [zeros(1,obj.M-1)+1,zeros(1,obj.D-obj.M+1)+2];
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
            K = ceil(obj.theta*(obj.D-obj.M+1));
            t = floor(obj.FE/obj.N/obj.taut)/obj.nt;
            g = sum(g1(X(:,obj.M:obj.M+K-1),0.25+0.75*sin(0.5*pi*t)),2);
            for i = 1 : ceil((obj.D-obj.M-K+1)/10)
                temp = 50 - sum(g1(X(:,obj.M+K+(i-1)*10:min(obj.M+K+i*10-1,end)),0),2);
                g(temp<50) = g(temp<50) + temp(temp<50);
            end
            PopObj = repmat(1+g/(obj.D-obj.M+1),1,obj.M).*fliplr(cumprod([ones(size(X,1),1),X(:,1:obj.M-1)],2)).*[ones(size(X,1),1),1-X(:,obj.M-1:-1:1)];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,1,10)';
                R = {a*a',a*(1-a'),(1-a)*ones(size(a'))};
            else
                R = [];
            end
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
            G      = 0.25+0.75*sin(0.5*pi*t);
            G      = round(G*1e6)/1e6;
            change = [0;find(G(1:end-1)~=G(2:end));length(G)];
            Scores = zeros(1,length(change)-1);
            for i = 1 : length(change)-1
                subPop    = Population(change(i)+1:change(i+1));
                Scores(i) = feval(metName,subPop,obj.optimum);
            end
            score = mean(Scores);
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
            G      = 0.25+0.75*sin(0.5*pi*t);
            G      = round(G*1e6)/1e6;
            change = [0;find(G(1:end-1)~=G(2:end));length(G)];
            tempStream = RandStream('mlfg6331_64','Seed',2);
            for i = 1 : length(change)-1
                color = rand(tempStream,1,3);
                Draw(Population(change(i)+1:change(i+1)).objs+(i-1)*0.1,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                Draw(obj.PF+(i-1)*0.1,'-','LineWidth',1,'Color',color);
            end
        end        
    end
end

function g = g1(x,t)
    g = (x-t).^2;
end
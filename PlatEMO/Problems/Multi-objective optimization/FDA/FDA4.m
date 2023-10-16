classdef FDA4 < PROBLEM
% <multi/many> <real> <large/none> <dynamic>
% Benchmark dynamic MOP proposed by Farina, Deb, and Amato
% taut --- 10 --- Number of generations for static optimization
% nt   --- 10 --- Number of distinct steps

%------------------------------- Reference --------------------------------
% M. Farina, K. Deb, and P. Amato, Dynamic multiobjective optimization
% problems: Test cases, approximations, and applications, IEEE Transactions
% on Evolutionary Computation, 2004, 8(5): 425-442.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        taut;       % Number of generations for static optimization
        nt;         % Number of distinct steps
        Optimums;   % Point sets on all Pareto fronts
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.taut,obj.nt] = obj.ParameterSet(10,10);
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+9; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
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
        function PopObj = CalObj(obj,PopDec)
            t = floor(obj.FE/obj.N/obj.taut)/obj.nt;
            G = abs(sin(0.5*pi*t));
            g = sum((PopDec(:,obj.M:end)-G).^2,2);
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
            G      = abs(sin(0.5.*pi.*t));
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
            G      = abs(sin(0.5.*pi.*t));
            G      = round(G*1e6)/1e6;
            change = [0;find(G(1:end-1)~=G(2:end));length(G)];
            tempStream = RandStream('mlfg6331_64','Seed',2);
            if obj.M == 2
                for i = 1 : length(change)-1
                    color = rand(tempStream,1,3);
                    Draw(Population(change(i)+1:change(i+1)).objs+(i-1)*0.1,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                    Draw(obj.PF+(i-1)*0.1,'-','LineWidth',1,'Color',color);
                end
            elseif obj.M == 3
                for i = 1 : length(change)-1
                    color = rand(tempStream,1,3);
                    ax = Draw(Population(change(i)+1:change(i+1)).objs+(i-1)*0.1,'o','MarkerSize',6,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
                    surf(ax,obj.PF{1}+(i-1)*0.1,obj.PF{2}+(i-1)*0.1,obj.PF{3}+(i-1)*0.1,'EdgeColor',color,'FaceColor','none');
                end
            else
                for i = 1 : length(change)-1
                    Draw(Population(change(i)+1:change(i+1)).objs,'-','Color',rand(tempStream,1,3),'LineWidth',2);
                end
            end
        end
    end
end
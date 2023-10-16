classdef FDA2 < PROBLEM
% <multi> <real> <large/none> <dynamic>
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
            obj.M = 2;
            if isempty(obj.D); obj.D = 10; end
            obj.D        = ceil((obj.D-1)/2)*2 + 1;
            obj.lower    = [0,-ones(1,obj.D-1)];
            obj.upper    = [1, ones(1,obj.D-1)];
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
            PopObj(:,1) = PopDec(:,1); 
            t = floor(obj.FE/obj.N/obj.taut)/obj.nt;
            g = 1 + sum(PopDec(:,2:(end-1)/2).^2,2);
            H = 0.75 + 0.7*sin(0.5.*pi.*t);
            % Note: The original definition of h is questionable
            h = 1 - (PopObj(:,1)./g).^(H+sum((PopDec(:,(end-1)/2+1:end)-H).^2,2));
            PopObj(:,2) = g.*h;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            % Generate point sets on all Pareto fronts
            t = floor((0:obj.maxFE)/obj.N/obj.taut)/obj.nt;
            H = 0.75 + 0.7*sin(0.5.*pi.*t);
            H = unique(round(H*1e6)/1e6);
            x = linspace(0,1,N)';
            obj.Optimums = {};
            for i = 1 : length(H)
                obj.Optimums(i,:) = {H(i),[x,1-x.^(H(i)+max(0,H(i)-1).^2*(obj.D-1)/2)]};
            end
            % Combine all point sets
            R = cat(1,obj.Optimums{:,2});
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
            H      = 0.75 + 0.7*sin(0.5.*pi.*t);
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
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
            H      = 0.75 + 0.7*sin(0.5.*pi.*t);
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
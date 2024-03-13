classdef MMF2 < PROBLEM
% <multi> <real> <multimodal>
% Multi-modal multi-objective test function

%------------------------------- Reference --------------------------------
% C. Yue, B. Qu, and J. Liang, A multi-objective particle swarm optimizer
% using ring topology for solving multimodal multiobjective Problems, IEEE
% Transactions on Evolutionary Computation, 2018, 22(5): 805-817.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        POS;    % Pareto optimal set for IGDX calculation
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [0,0];
            obj.upper    = [1,2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            temp = X(:,2) <= 1;
            y    = zeros(size(X,1),1);
            y(temp)  = X(temp,2) - sqrt(X(temp,1));
            y(~temp) = X(~temp,2) - 1 - sqrt(X(~temp,1));
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = 1 - sqrt(X(:,1)) + 2.*((4.*y.^2)-2*cos(20.*y.*pi/sqrt(2))+2);
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj,N)
            % Generate points in Pareto optimal set
            obj.POS(:,1) = linspace(0,1,N/2)';
            obj.POS(:,2) = sqrt(obj.POS(:,1));
            obj.POS = [obj.POS;obj.POS(:,1),obj.POS(:,2)+1];
            % Generate points on Pareto front
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - sqrt(R(:,1));
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R(:,1) = linspace(0,1,100)';
            R(:,2) = 1 - sqrt(R(:,1));
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
            temp   = PopDec(:,2)<=1;
            Draw(Population(temp).objs,'o','MarkerSize',6,'Marker','o','Markerfacecolor',[1 .5 .5],'Markeredgecolor',[1 .2 .2],{'\it f\rm_1','\it f\rm_2',[]});
            Draw(Population(~temp).objs+0.1,'o','MarkerSize',6,'Marker','o','Markerfacecolor',[.5 .5 1],'Markeredgecolor',[.2 .2 1]);
            Draw(obj.PF,'-','LineWidth',1,'Color',[1 .2 .2]);
            Draw(obj.PF+0.1,'-','LineWidth',1,'Color',[.2 .2 1]);
        end
    end
end
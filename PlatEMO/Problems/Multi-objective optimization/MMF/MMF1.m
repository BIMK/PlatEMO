classdef MMF1 < PROBLEM
% <multi> <real> <multimodal>
% Multi-modal multi-objective test function

%------------------------------- Reference --------------------------------
% C. Yue, B. Qu, and J. Liang, A multi-objective particle swarm optimizer
% using ring topology for solving multimodal multiobjective Problems, IEEE
% Transactions on Evolutionary Computation, 2018, 22(5): 805-817.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [1,-1];
            obj.upper    = [3,1];
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            PopObj(:,1) = abs(X(:,1)-2);
            PopObj(:,2) = 1-sqrt(PopObj(:,1))+2*(X(:,2)-sin(6*pi*PopObj(:,1)+pi)).^2; 
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(1,3,N)';
            R(:,2) = sin(6*pi*abs(R(:,1)-2)+pi);
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            PopDec = Population.decs;
            temp   = PopDec(:,1)<=2;
            Draw(Population(temp).objs,'o','MarkerSize',6,'Marker','o','Markerfacecolor',[1 .5 .5],'Markeredgecolor',[1 .2 .2],{'\it f\rm_1','\it f\rm_2',[]});
            Draw(Population(~temp).objs+0.1,'o','MarkerSize',6,'Marker','o','Markerfacecolor',[.5 .5 1],'Markeredgecolor',[.2 .2 1]);
            L  = [0:0.01:1;1-sqrt(0:0.01:1)]';
            Draw(L,'-','LineWidth',1,'Color',[1 .2 .2]);
            Draw(L+0.1,'-','LineWidth',1,'Color',[.2 .2 1]);
        end
    end
end
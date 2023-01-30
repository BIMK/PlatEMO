classdef MPDMP < PROBLEM
% <multi/many> <real>
% The multi-point distance minimization problem
% lower --- -100 --- Lower bound of decision variables
% upper ---  100 --- Upper bound of decision variables

%------------------------------- Reference --------------------------------
% M. Koppen and K. Yoshida, Substitute distance assignments in NSGA-II for
% handling many-objective optimization problems, Proceedings of the
% International Conference on Evolutionary Multi-Criterion Optimization,
% 2007, 727-741.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Points; % Vertexes
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [lower,upper] = obj.ParameterSet(-100,100);
            if isempty(obj.M); obj.M = 10; end
            obj.M        = max(obj.M,3);
            obj.D        = 2;
            obj.lower    = zeros(1,2) + lower;
            obj.upper    = zeros(1,2) + upper;
            obj.encoding = ones(1,obj.D);
            % Generate vertexes
            if mod(obj.M,2) == 0
                Angle = (2.*(1:obj.M)-3).*pi./obj.M;
            else
                Angle = (2.*(1:obj.M)-2).*pi./obj.M;
            end
            obj.Points = [sin(Angle)',cos(Angle)'];
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = pdist2(PopDec,obj.Points);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            [X,Y] = ndgrid(linspace(-1,1,ceil(sqrt(N))));
            ND    = inpolygon(X(:),Y(:),obj.Points(:,1),obj.Points(:,2));
            R     = pdist2([X(ND),Y(ND)],obj.Points);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 3
                [X,Y]    = ndgrid(linspace(-1,1,40));
                R        = pdist2([X(:),Y(:)],obj.Points);
                ND       = inpolygon(X(:),Y(:),obj.Points(:,1),obj.Points(:,2));
                R(~ND,:) = nan;
                R = {reshape(R(:,1),size(X)),reshape(R(:,2),size(X)),reshape(R(:,3),size(X))};
            else
                R = [];
            end
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            Draw(obj.Points([1:end,1],:),'-k','LineWidth',1.5,{'\it x\rm_1','\it x\rm_2',[]});
            Draw(obj.Points,'o','MarkerSize',6,'Marker','o','Markerfacecolor',[1 1 1],'Markeredgecolor',[.4 .4 .4]);
            Draw(Population.decs);
        end
    end
end
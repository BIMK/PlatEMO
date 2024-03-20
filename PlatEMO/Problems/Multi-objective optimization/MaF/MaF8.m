classdef MaF8 < PROBLEM
% <multi/many> <real>
% MP-DMP

%------------------------------- Reference --------------------------------
% R. Cheng, M. Li, Y. Tian, X. Zhang, S. Yang, Y. Jin, and X. Yao, A
% benchmark test suite for evolutionary many-objective optimization,
% Complex & Intelligent Systems, 2017, 3(1): 67-81.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
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
            if isempty(obj.M); obj.M = 10; end
            obj.M        = max(obj.M,3);
            obj.D        = 2;
            obj.lower    = [-10000,-10000];
            obj.upper    = [10000,10000];
            obj.encoding = ones(1,obj.D);
            % Generate vertexes
            obj.Points  = [];
            [thera,rho] = cart2pol(0,1);
            [obj.Points(:,1),obj.Points(:,2)] = pol2cart(thera-(1:obj.M)*2*pi/obj.M,rho);
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
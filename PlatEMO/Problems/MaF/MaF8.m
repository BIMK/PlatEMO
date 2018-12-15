classdef MaF8 < PROBLEM
% <problem> <MaF>
% MP-DMP

%------------------------------- Reference --------------------------------
% R. Cheng, M. Li, Y. Tian, X. Zhang, S. Yang, Y. Jin, and X. Yao, A
% benchmark test suite for evolutionary many-objective optimization,
% Complex & Intelligent Systems, 2017, 3(1): 67-81.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
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
        %% Initialization
        function obj = MaF8()
            % Parameter setting
            if isempty(obj.Global.M)
                obj.Global.M = 10;
            end
            obj.Global.D        = 2;
            obj.Global.lower    = [-10000,-10000];
            obj.Global.upper    = [10000,10000];
            obj.Global.encoding = 'real';
            % Generate vertexes
            obj.Points  = [];
            [thera,rho] = cart2pol(0,1);
            [obj.Points(:,1),obj.Points(:,2)] = pol2cart(thera-(1:obj.Global.M)*2*pi/obj.Global.M,rho);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = pdist2(PopDec,obj.Points);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            [X,Y] = ndgrid(linspace(-1,1,ceil(sqrt(N))));
            ND    = inpolygon(X(:),Y(:),obj.Points(:,1),obj.Points(:,2));
            P     = pdist2([X(ND),Y(ND)],obj.Points);
        end
        %% Draw special figure
        function Draw(obj,PopDec)
            cla; Draw(PopDec);
            plot(obj.Points([1:end,1],1),obj.Points([1:end,1],2),'-k','LineWidth',1.5);
            plot(obj.Points(:,1),obj.Points(:,2),'ok','MarkerSize',6,'Marker','o','Markerfacecolor',[1 1 1],'Markeredgecolor',[.4 .4 .4]);
            xlabel('\itx\rm_1'); ylabel('\itx\rm_2');
        end            
    end
end
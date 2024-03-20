classdef MaF7 < PROBLEM
% <multi/many> <real> <large/none>
% DTLZ7

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

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M + 19; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),obj.M);
            g      = 1+9*mean(PopDec(:,obj.M:end),2);
            PopObj(:,1:obj.M-1) = PopDec(:,1:obj.M-1);
            PopObj(:,obj.M)     = (1+g).*(obj.M-sum(PopObj(:,1:obj.M-1)./(1+repmat(g,1,obj.M-1)).*(1+sin(3*pi.*PopObj(:,1:obj.M-1))),2));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            interval     = [0,0.251412,0.631627,0.859401];
            median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
            X            = UniformPoint(N,obj.M-1,'grid');
            X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
            X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
            R            = [X,2*(obj.M-sum(X/2.*(1+sin(3*pi.*X)),2))];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                x      = linspace(0,1,100)';
                y      = 2*(2-x/2.*(1+sin(3*pi*x)));
                nd     = NDSort([x,y],1)==1;
                x(~nd) = nan;
                R      = [x,y];
            elseif obj.M == 3
                [x,y]  = meshgrid(linspace(0,1,20));
                z      = 2*(3-x/2.*(1+sin(3*pi*x))-y/2.*(1+sin(3*pi*y)));
                nd     = reshape(NDSort([x(:),y(:),z(:)],1)==1,size(z));
                z(~nd) = nan;
                R      = {x,y,z};
            else
                R = [];
            end
        end
    end
end
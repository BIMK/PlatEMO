classdef MinusDTLZ7 < PROBLEM
% <2017> <multi/many> <real> <large/none> <expensive/none>
% Minus-DTLZ7

%------------------------------- Reference --------------------------------
% H. Ishibuchi, Y. Setoguchi, H. Masuda, and Y. Nojima. Performance of
% decomposition-based many-objective algorithms strongly depends on Pareto
% front shapes. IEEE Transactions on Evolutionary Computation, 2017, 21(2):
% 169-190.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yuma Horaguchi and Ryudai Kato (email: horaguchi-yuma-hv@ynu.jp)

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+19; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),obj.M);
            g      = 1+9*mean(PopDec(:,obj.M:end),2);
            PopObj(:,1:obj.M-1) = -1*PopDec(:,1:obj.M-1);
            PopObj(:,obj.M)     = -1*(1+g).*(obj.M-sum(PopObj(:,1:obj.M-1)./(1+repmat(g,1,obj.M-1)).*(1+sin(3*pi.*PopObj(:,1:obj.M-1))),2));
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            N_split      = ceil(N ^ (1 / (obj.M - 1)));
            area_min     = 0.541091;
            area_max     = 0.581368;
            out_of_area  = (area_max - area_min) * ((N_split - 1) / (N_split - 2)) + area_min;
            X            = UniformPoint(N,obj.M-1,'grid');
            X            = X .* (out_of_area - area_min) + area_min;
            X(X > area_max + 1e-6) = 1;
            g            = 10;
            R            = [X,(1+g)*(obj.M+sum(X/(1+g).*(1-sin(3*pi.*X)),2))];
            R            = R*(-1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                g      = 10;
                x      = linspace(0,1,201)';
                y      = (1+g)*(2+x/(1+g).*(1-sin(3*pi*x)));
                nd     = NDSort([-x,-y],1)==1;
                x(~nd) = nan;
                R      = [-x,-y];
            elseif obj.M == 3
                g      = 10;
                [x,y]  = meshgrid(linspace(0,1,51));
                z      = (1+g)*(3+x/(1+g).*(1-sin(3*pi*x))+y/(1+g).*(1-sin(3*pi*y)));
                nd     = reshape(NDSort([-x(:),-y(:),-z(:)],1)==1,size(z));
                z(~nd) = nan;
                R      = {-x,-y,-z};
            else
                R = [];
            end
        end
    end
end
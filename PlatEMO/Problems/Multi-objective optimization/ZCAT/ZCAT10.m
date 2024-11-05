classdef ZCAT10 < PROBLEM
% <2023> <multi/many> <real> <large/none> <expensive/none>
% Benchmark MOP proposed by Zapotecas, Coello, Aguirre, and Tanaka
% complicatedPS --- 0 --- Whether complicated Pareto set is considered
% bias          --- 0 --- Whether bias is considered
% imbalance     --- 0 --- Whether imbalance is considered

%------------------------------- Reference --------------------------------
% S. Zapotecas-Martinez, C. A. Coello Coello, H. E. Aguirre, and K. Tanaka,
% Challenging test problems for multi- and many-objective optimization,
% Swarm and Evolutionary Computation, 2023, 81: 101350.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        complicatedPS = 0;
        bias          = 0;
        imbalance     = 0;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.complicatedPS,obj.bias,obj.imbalance] = obj.ParameterSet(0,0,0);
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = 2*obj.M; end
            obj.D        = max(obj.D,2*obj.M);
            obj.lower    = -(1:obj.D)/2;
            obj.upper    = (1:obj.D)/2;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            y = (PopDec-repmat(obj.lower,size(PopDec,1),1))./repmat(obj.upper-obj.lower,size(PopDec,1),1);
            A = repmat((1:obj.M).^2,size(y,1),1).*[y(:,1:obj.M-1),(1/0.02-1./(mean(1-y(:,1:obj.M-1),2)+0.02))./(1/0.02-1/1.02)];
            if obj.complicatedPS
                func = @g7;
            else
                func = @g0;
            end
            for j = obj.M : obj.D
                z(:,j) = y(:,j) - func(y(:,1:obj.M-1),2*pi*(j-obj.M)/obj.D);
            end
            if obj.bias
                w = abs(z).^0.05;
            else
                w = abs(z);
            end
            for i = 1 : obj.M
                if obj.imbalance
                    B(:,i) = i^2*10/(2*exp(1)-2)*(exp(max(sqrt(abs(w(:,obj.M-1+i:obj.M:end))),[],2))-exp(mean(1/2*(cos(9*pi*w(:,obj.M-1+i:obj.M:end))+1),2))-1+exp(1));
                else
                    B(:,i) = i^2*10*mean(w(:,obj.M-1+i:obj.M:end).^2,2);
                end
            end
            PopObj = A + B;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            x = UniformPoint(N,obj.M-1,'grid');
            R = [x,(1/0.02-1./(mean(1-x,2)+0.02))./(1/0.02-1/1.02)];
            R = repmat((1:obj.M).^2,size(R,1),1).*R;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                [x,y] = meshgrid(linspace(0,1,10));
                R     = {x,4*y,9*(1/0.02-1./(1-x/2-y/2+0.02))./(1/0.02-1/1.02)};
            else
                R = [];
            end
        end
    end
end

function g = g0(y,~)
    g = zeros(size(y,1),1) + 0.2210;
end

function g = g7(y,theta)
    g = (mean(y,2)+exp(sin(7*pi*mean(y,2)-pi/2+theta))-exp(-1))./(1+exp(1)-exp(-1));
end
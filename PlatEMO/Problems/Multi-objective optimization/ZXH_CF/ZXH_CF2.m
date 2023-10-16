classdef ZXH_CF2 < PROBLEM
% <multi/many> <real> <large/none> <constrained>
% Constrained benchmark MOP proposed by Zhou, Xiang, and He

%------------------------------- Reference --------------------------------
% Y. Zhou, Y. Xiang, and X. He, Constrained multiobjective optimization:
% Test problem construction and performance evaluations, IEEE Transactions
% on Evolutionary Computation, 2021, 25(1): 172-186.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
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
            if isempty(obj.D); obj.D = obj.M+10;  end
            obj.lower    = zeros(1,obj.D) + 1e-10;
            obj.upper    = ones(1,obj.D)  - 1e-10;
            obj.encoding = ones(1,obj.D); 
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            PopDec = varargin{1};
            OptX   = 0.2;               
            [N,D]  = size(PopDec);
            M      = obj.M;
            % Step 1: Compute cumsum 
            Sx = cumsum(PopDec(:,1:M).^2,2,'reverse'); 
            % Step 2: Compute theta
            THETA = 2/pi*atan(sqrt(Sx(:,2:end))./PopDec(:,1:M-1));
            % Step 3: Calculate Rosenbrock function
            h = sum(100*((PopDec(:,M+1:end-1)-OptX).^2-(PopDec(:,M+2:end)-OptX)).^2+(PopDec(:,M+1:end-1)-OptX).^2,2);
            % Step 4: Compute T_
            T = (1 - Sx(:,1)).^2 + h;
            % Step 5: Objectives (concave)
            G      = [ones(N,1) cumprod(sin(pi/2*THETA),2)] .* [cos(pi/2*THETA) ones(N,1)];
            PopObj = G .* repmat((1+T),1,M) ;          
            % Step 6: Constraints
            PopCon(:,1) = Sx(:,1) + h - 1;
            PopCon(:,2) = -(Sx(:,1) + h - 1/4);
            Population  = SOLUTION(varargin{1},PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
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
                R = {cos(a)*cos(a'),cos(a)*sin(a'),sin(a)*ones(size(a'))};
            else
                R = [];
            end
        end
    end
end